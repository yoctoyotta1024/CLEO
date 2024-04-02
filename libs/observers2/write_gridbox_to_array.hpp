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
 * Last Modified: Tuesday 2nd April 2024
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
#include "./xarray_for_gridbox_data.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/**
 * @brief Concept for WriteGridboxToArray is all types that have functions for creating a functor
 * to collect data from a gridbox (in a parallel for loop) and then write the collected data for all
 * gridboxes to an array in a dataset.
 *
 * @tparam WG2A The type that satisfies the WriteGridboxToArray concept.
 */
template <typename WG2A, typename Store>
concept WriteGridboxToArray =
    requires(WG2A wg2a, const Dataset<Store> &ds, const viewd_constgbx d_gbxs) {
      { wg2a.get_functor(d_gbxs) };
      { wg2a.write_to_array(ds) } -> std::same_as<void>;
      { wg2a.write_arrayshape(ds) } -> std::same_as<void>;
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

    /* Functor operator to perform copy of each element in parallel in Kokkos Range Policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const {
      a_functor(ii);
      b_functor(ii);
    }

    /* Functor operator to perform copy of each element in parallel in Kokkos Team Policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {
      a_functor(team_member);
      b_functor(team_member);
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

  void write_to_array(const Dataset<Store> &dataset) const {
    a.write_to_array(dataset);
    b.write_to_array(dataset);
  }

  void write_arrayshape(const Dataset<Store> &dataset) const {
    a.write_arrayshape(dataset);
    b.write_arrayshape(dataset);
  }
};

/* struct satifying WriteGridboxToArray and does nothing */
template <typename Store>
struct NullWriteGridboxToArray {
 public:
  struct Functor {
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const {}
  };

  Functor get_functor(const viewd_constgbx d_gbxs) const { return Functor{}; }

  void write_to_array(const Dataset<Store> &dataset) const {}

  void write_arrayshape(const Dataset<Store> &dataset) const {}
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
struct CombineWG2A {
  template <WriteGridboxToArray<Store> WriteGbx1, WriteGridboxToArray<Store> WriteGbx2>
  auto operator()(const WriteGbx1 a, const WriteGbx2 b) const {
    return CombinedWriteGridboxToArray<Store, WriteGbx1, WriteGbx2>(a, b);
  }
};

#endif  // LIBS_OBSERVERS2_WRITE_GRIDBOX_TO_ARRAY_HPP_
