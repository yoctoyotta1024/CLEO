/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: write_gridbox_to_array.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 29th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Concept and structs to write data collected from each Gridbox in parallel
 * to an array in a dataset
 */

#ifndef LIBS_OBSERVERS2_WRITE_GRIDBOX_TO_ARRAY_HPP_
#define LIBS_OBSERVERS2_WRITE_GRIDBOX_TO_ARRAY_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

/**
 * @brief Concept for WriteGridboxToArray is all types that have functions for creating a functor
 * to collect data from a gridbox (in a parallel for loop) and then write the collected data for all
 * gridboxes to an array in a dataset.
 *
 * @tparam GDW The type that satisfies the WriteGridboxToArray concept.
 */
template <typename GDW, typename Store>
concept WriteGridboxToArray = requires(GDW gdw, Dataset<Store> &ds, const viewd_constgbx d_gbxs) {
  { gdw.get_functor(d_gbxs) };
  { gdw.write_to_array(ds) } -> std::same_as<void>;
  { gdw.write_arrayshape(ds) } -> std::same_as<void>;
};

/**
 * @brief Combined gridbox data writer struct combines two structs that write gridbox data to an
 * array into one struct that does the actions of both.
 *
 * @tparam WriteGbx1 The type of the first gridbox data writer.
 * @tparam WriteGbx2 The type of the second gridbox data writer.
 */
template <typename Store, WriteGridboxToArray<Store> WriteGbx1,
          WriteGridboxToArray<Store> WriteGbx2>
struct CombinedWriteGridboxToArray {
 private:
  WriteGbx1 a; /**< The first instance of type of WriteGridboxToArray. */
  WriteGbx2 b; /**< The second instance of type of WriteGridboxToArray. */

 public:
  struct Functor {
    WriteGbx1::Functor a_functor;
    WriteGbx2::Functor b_functor;

    explicit Functor(const WriteGbx1 a, const WriteGbx2 b, const viewd_constgbx d_gbxs)
        : a_functor(a.get_functor(d_gbxs)), b_functor(b.get_functor(d_gbxs)) {}

    // Functor operator to perform copy of each element in parallel
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const {
      a_functor(ii);
      b_functor(ii);
    }
  };

  /**
   * @brief Constructs a CombinedWriteGridboxToArray object.
   *
   * @param a The first gridbox data writer.
   * @param b The second gridbox data writer.
   */
  CombinedWriteGridboxToArray(const WriteGbx1 a, const WriteGbx2 b) : a(a), b(b) {}

  Functor get_functor(const viewd_constgbx d_gbxs) const { return Functor(a, b, d_gbxs); }

  void write_to_array(Dataset<Store> &dataset) const {
    a.write_to_array(dataset);
    b.write_to_array(dataset);
  }

  void write_arrayshape(Dataset<Store> &dataset) const {
    a.write_arrayshape(dataset);
    b.write_arrayshape(dataset);
  }
};

/**
 * @brief Operator for combining two gridbox data writers which write to a FSStore.
 *
 * This operator combines two gridbox data writers into one using the
 * CombinedWriteGridboxToArray struct with Store = FSStore.
 *
 * @param a The first gridbox data writer.
 * @param b The second gridbox data writer.
 * @return The combined gridbox data writer.
 */
template <typename Store>
struct CombineGDW {
  template <WriteGridboxToArray<Store> WriteGbx1, WriteGridboxToArray<Store> WriteGbx2>
  auto operator()(const WriteGbx1 a, const WriteGbx2 b) const {
    return CombinedWriteGridboxToArray<Store, WriteGbx1, WriteGbx2>(a, b);
  }
};

// struct satifying WriteGridboxToArray and does nothing
template <typename Store>
struct NullGbxWriter {
 public:
  struct Functor {
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const {}
  };

  Functor get_functor(const viewd_constgbx d_gbxs) const { return Functor{}; }

  void write_to_array(Dataset<Store> &dataset) const {}

  void write_arrayshape(Dataset<Store> &dataset) const {}
};

/* struct holding an array in a dataset as well a view and its mirror view
which can be useful when collecting data for 1 variable from 'ngbxs' gridboxes
(in parallel) to then write the the array */
template <typename Store, typename T>
struct XarrayForGenericGbxWriter {
  XarrayZarrArray<Store, T> xzarr;                         // array in a dataset
  using viewh_data = Buffer<T>::viewh_buffer;              // type of view for h_data
  using mirrorviewd_data = Buffer<T>::mirrorviewd_buffer;  // mirror view type for d_data
  viewh_data h_data;        // view on host for value of 1 variable from every gridbox
  mirrorviewd_data d_data;  // mirror view of h_data on device

  // Constructor to initialize views and pointer to array in dataset
  XarrayForGenericGbxWriter(Dataset<Store> &dataset, const std::string_view name,
                            const std::string_view units, const std::string_view dtype,
                            const double scale_factor, const size_t maxchunk, const size_t ngbxs)
      : xzarr(dataset.template create_array<T>(name, units, dtype, scale_factor,
                                               good2Dchunkshape(maxchunk, ngbxs),
                                               {"time", "gbxindex"})),
        h_data("h_data", ngbxs),
        d_data(Kokkos::create_mirror_view(ExecSpace(), h_data)) {}

  // copy data from device view directly to host and then write to array in dataset
  void write_to_array(Dataset<Store> &dataset) {
    Kokkos::deep_copy(h_data, d_data);
    dataset.write_to_array(xzarr, h_data);
  }

  // call function to write shape of array according to dataset
  void write_arrayshape(Dataset<Store> &dataset) { dataset.write_arrayshape(xzarr); }
}

// template WriteGridboxToArray to write one variable from each gridbox to an array in a dataset
template <typename Store, typename T, typename FunctorFunc>
class GenericGbxWriter {
 private:
  std::shared_ptr<XarrayForGenericGbxWriter<Store, T>> xzarr_ptr;
  FunctorFunc ffunc;

 public:
  struct Functor {
    using mirrorviewd_data = XarrayForGenericGbxWriter<Store, T>::mirrorviewd_data;
    FunctorFunc ffunc;
    viewd_constgbx d_gbxs;    // view of gridboxes on device
    mirrorviewd_data d_data;  // mirror view for data on device

    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, mirrorviewd_data d_data)
        : ffunc(ffunc), d_gbxs(d_gbxs), d_data(d_data) {}

    // Functor operator to perform copy of 1 variable in each gridbox to d_data in parallel
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const { ffunc(ii, d_gbxs, d_data); }
  };

  // Constructor to initialize views and pointer to array in dataset
  GenericGbxWriter(std::shared_ptr<XarrayForGenericGbxWriter<Store, T>> xzarr_ptr,
                   FunctorFunc ffunc)
      : xzarr_ptr(xzarr_ptr), ffunc(ffunc) {}

  // return functor for getting 1 variable from every gridbox in parallel
  Functor get_functor(const viewd_constgbx d_gbxs) const {
    assert((d_gbxs.extent(0) == xzarr_ptr->d_data.extent(0)) &&
           "d_data view must be size of the number of gridboxes");
    return Functor(ffunc, d_gbxs, xzarr_ptr->d_data);
  }

  // copy data from device view directly to host and then write to array in dataset
  void write_to_array(Dataset<Store> &dataset) const { xzarr_ptr->write_to_array(dataset); }

  // call function to write shape of array according to dataset
  void write_arrayshape(Dataset<Store> &dataset) const { xzarr_ptr->write_arrayshape(dataset); }
};

#endif  // LIBS_OBSERVERS2_WRITE_GRIDBOX_TO_ARRAY_HPP_
