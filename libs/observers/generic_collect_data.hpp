/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: generic_collect_data.hpp
 * Project: observers
 * Created Date: Thursday 4th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Very generic struct satisyfing the CollectDataForDataset concept to collect data for a
 * variable from gridboxes and/or superdroplets and write it to an xarray in a datatset.
 */

#ifndef LIBS_OBSERVERS_GENERIC_COLLECT_DATA_HPP_
#define LIBS_OBSERVERS_GENERIC_COLLECT_DATA_HPP_

#include <Kokkos_Core.hpp>
#include <memory>

#include "../kokkosaliases.hpp"
#include "zarr/buffer.hpp"
#include "zarr/xarray_zarr_array.hpp"

/**
 * @brief Struct to 1) manage collecting data into a view in host memory by copying data from the
 * device execution space and 2) manage how to write this data to an Xarray for a variable in a
 * dataset.
 *
 * This struct manages an Xarray, a view in host memory and a mirror view in device memory in order
 * to collect data for a variable from the execution space and write it to an Xarray in a dataset.
 *
 * @tparam Store The type of the data store in the Xarray.
 * @tparam T The type of the data in the Xarray.
 */
template <typename Store, typename T>
struct XarrayAndViews {
  using viewh_data = Buffer<T>::viewh_buffer;             /**< type of view for h_data */
  using mirrorviewd_data = Buffer<T>::mirrorviewd_buffer; /**< mirror view type for d_data */
  XarrayZarrArray<Store, T> xzarr; /**< Xarray with Zarr backend to write h_data to */
  viewh_data h_data;               /**< view on host used to collect some data for the Xarray */
  mirrorviewd_data d_data;         /**< mirror view of h_data on device */

  /**
   * @brief Constructs a new XarrayAndViews object.
   *
   * @param xzarr The Xarray with Zarr backend object.
   * @param dataview_size The size of the views for collecting data.
   */
  XarrayAndViews(const XarrayZarrArray<Store, T> xzarr, const size_t dataview_size)
      : xzarr(xzarr),
        h_data("h_data", dataview_size),
        d_data(Kokkos::create_mirror_view(ExecSpace(), h_data)) {}
};

/* generic struct satisyfing the CollectDataForDataset concept to collect data for a
 * variable and write it to an xarray in a datatset. */

/**
 * @brief Generic class satisyfing the CollectDataForDataset concept to collect data for a variable
 * and write it to an Xarray in a dataset.
 *
 * This class provides a functor to collect data into a view in device memory for a single variable
 * from superdroplets and/or gridboxes (to be used in a parallel range policy loop over gridboxes
 * and/or superdroplets). It also provides functions to then write that collected data to an Xarray
 * in a dataset.
 *
 * @tparam Store The type of the data store of the dataset.
 * @tparam T The type of the data in the Xarray.
 * @tparam FunctorFunc Type to act as a functor in a Kokkos parallel range policy loop for
 * collecting data for a variable from gridbxoes and/or superdroplets
 */
template <typename Store, typename T, typename FunctorFunc>
class GenericCollectData {
 private:
  FunctorFunc ffunc; /**< functor to collect data into a view during a parallel range policy loop */
  std::shared_ptr<XarrayAndViews<Store, T>> ptr; /**< pointer to xarray and views to collect data */

 public:
  /**
   * @brief Generic wrapper to use FunctorFunc type to collect data into a view in device memory
   * during a Kokkos::parallel_for loop with a range policy.
   */
  struct Functor {
    using mirrorviewd_data = XarrayAndViews<Store, T>::mirrorviewd_data;
    FunctorFunc ffunc;             /**< functor to collect data into d_data during parallel loop */
    viewd_constgbx d_gbxs;         /**< view of gridboxes on device */
    subviewd_constsupers d_supers; /**< view of superdroplets on device */
    mirrorviewd_data d_data;       /**< mirror view on device for data to collect */

    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
            mirrorviewd_data d_data)
        : ffunc(ffunc), d_gbxs(d_gbxs), d_supers(d_supers), d_data(d_data) {}

    /**
     * @brief Adapter from signature of Kokkos::parallel_for with a range policy to call to
     * FunctorFunc type for collecting data into d_data from gridboxes and/or superdroplets.
     *
     * @param nn The index of the data element.
     */
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t nn) const { ffunc(nn, d_gbxs, d_supers, d_data); }
  };

  /* Constructor to initialize GenericCollectData given functor function-like object,
  an xarray in a dataset and the size of the data view  */
  /**
   * @brief Constructs a new GenericCollectData object.
   *
   * The dataview_size should match the number of elements to collect when ffunc is called during a
   * Kokkos::parallel_for loop which uses a range policy over gridboxes and/or superdroplets.
   *
   * @param ffunc Function-like object for the functor to collect data.
   * @param xzarr The Xarray object in a dataset to write for a variable data to.
   * @param dataview_size The size of the view to collect data (number of elements).
   */
  GenericCollectData(const FunctorFunc ffunc, const XarrayZarrArray<Store, T> xzarr,
                     const size_t dataview_size)
      : ffunc(ffunc), ptr(std::make_shared<XarrayAndViews<Store, T>>(xzarr, dataview_size)) {}

  /**
   * @brief Returns the functor for collecting data.
   *
   * @param d_gbxs The view of gridboxes on device.
   * @param d_supers The view of superdroplets on device.
   * @return The functor object to use during a Kokkos:parallel_for range policy loop.
   */
  Functor get_functor(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    // assert(((ptr->d_data.extent(0) == d_gbxs.extent(0)) ||
    //         (ptr->d_data.extent(0) == d_supers.extent(0))) &&
    //        "d_data view should be size of the number of gridboxes or superdroplets");
    return Functor(ffunc, d_gbxs, d_supers, ptr->d_data);
  }

  /**
   * @brief Reallocates the views with a new size.
   *
   * The size of the view should match the number of elements to collect when ffunc is called during
   * a Kokkos::parallel_for loop which uses a range policy over gridboxes and/or superdroplets.
   *
   * @param size The new size of the views on device and host.
   */
  void reallocate_views(const size_t size) const {
    Kokkos::realloc(ptr->h_data, size);
    Kokkos::realloc(ptr->d_data, size);
  }

  /**
   * @brief Deep copies data for an array from the device view to the host and then writes it to an
   * array in the dataset.
   *
   * @param dataset The dataset to write data to.
   */
  template <typename Dataset>
  void write_to_arrays(const Dataset &dataset) const {
    Kokkos::deep_copy(ptr->h_data, ptr->d_data);
    dataset.write_to_array(ptr->xzarr, ptr->h_data);
  }

  /**
   * @brief Deep copies data for a ragged array from the device view to the host and then writes it
   * to a ragged array in the dataset.
   *
   * @param dataset The dataset to write data to.
   */
  template <typename Dataset>
  void write_to_ragged_arrays(const Dataset &dataset) const {
    Kokkos::deep_copy(ptr->h_data, ptr->d_data);
    dataset.write_to_ragged_array(ptr->xzarr, ptr->h_data);
  }

  /**
   * @brief Calls a function to write the shape of an array to the dataset.
   *
   * @param dataset The dataset to write array shape to.
   */
  template <typename Dataset>
  void write_arrayshapes(const Dataset &dataset) const {
    dataset.write_arrayshape(ptr->xzarr);
  }

  /**
   * @brief Calls a function to write the shape of a ragged array to the dataset.
   *
   * @param dataset The dataset to write array shape to.
   */
  template <typename Dataset>
  void write_ragged_arrayshapes(const Dataset &dataset) const {
    dataset.write_ragged_arrayshape(ptr->xzarr);
  }
};

#endif  // LIBS_OBSERVERS_GENERIC_COLLECT_DATA_HPP_
