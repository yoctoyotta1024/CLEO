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
 * @brief Combined gridbox data writer struct combines two gridbox data writers into one.
 *
 * @tparam GbxWriter1 The type of the first gridbox data writer.
 * @tparam GbxWriter2 The type of the second gridbox data writer.
 */
template <typename Store, WriteGridboxToArray<Store> GbxWriter1,
          WriteGridboxToArray<Store> GbxWriter2>
struct CombinedWriteGridboxToArray {
 private:
  GbxWriter1 a; /**< The first instance of type of WriteGridboxToArray. */
  GbxWriter2 b; /**< The second instance of type of WriteGridboxToArray. */

 public:
  struct Functor {
    GbxWriter1::Functor a_functor;
    GbxWriter2::Functor b_functor;

    explicit Functor(const GbxWriter1 a, const GbxWriter2 b, const viewd_constgbx d_gbxs)
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
  CombinedWriteGridboxToArray(const GbxWriter1 a, const GbxWriter2 b) : a(a), b(b) {}

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
  template <WriteGridboxToArray<Store> GbxWriter1, WriteGridboxToArray<Store> GbxWriter2>
  auto operator()(const GbxWriter1 a, const GbxWriter2 b) const {
    return CombinedWriteGridboxToArray<Store, GbxWriter1, GbxWriter2>(a, b);
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

// template WriteGridboxToArray to write one variable from each gridbox to an array in a dataset
template <typename Store, typename T, typename FunctorFunc>
class GenericGbxWriter {
 private:
  FunctorFunc ffunc;
  using viewh_data = Buffer<T>::viewh_buffer;              // type of view for h_data
  using mirrorviewd_data = Buffer<T>::mirrorviewd_buffer;  // mirror view type for d_data
  std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr;    // pointer to array in dataset
  viewh_data h_data;        // view on host for value of 1 variable from every gridbox
  mirrorviewd_data d_data;  // mirror view of h_data on device

 public:
  struct Functor {
    FunctorFunc ffunc;
    viewd_constgbx d_gbxs;    // view of gridboxes
    mirrorviewd_data d_data;  // mirror view on device for value of 1 variable from every gridbox

    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, mirrorviewd_data d_data)
        : ffunc(ffunc), d_gbxs(d_gbxs), d_data(d_data) {}

    // Functor operator to perform copy of 1 variable in each gridbox to d_data in parallel
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const { ffunc(ii, d_gbxs, d_data); }
  };

  // Constructor to initialize views and pointer to array in dataset
  GenericGbxWriter(Dataset<Store> &dataset, FunctorFunc ffunc,
                   std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr, const size_t ngbxs)
      : ffunc(ffunc),
        xzarr_ptr(xzarr_ptr),
        h_data("h_data", ngbxs),
        d_data(Kokkos::create_mirror_view(ExecSpace(), h_data)) {}

  // return functor for getting 1 variable from every gridbox in parallel
  Functor get_functor(const viewd_constgbx d_gbxs) const {
    assert((d_gbxs.extent(0) == d_data.extent(0)) &&
           "d_data view must be size of the number of gridboxes");
    return Functor(ffunc, d_gbxs, d_data);
  }

  // copy data from device view directly to host and then write to array in dataset
  void write_to_array(Dataset<Store> &dataset) const {
    Kokkos::deep_copy(h_data, d_data);
    dataset.write_to_array(xzarr_ptr, h_data);
  }

  // call function to write shape of array according to dataset
  void write_arrayshape(Dataset<Store> &dataset) const { dataset.write_arrayshape(xzarr_ptr); }
};

#endif  // LIBS_OBSERVERS2_WRITE_GRIDBOX_TO_ARRAY_HPP_
