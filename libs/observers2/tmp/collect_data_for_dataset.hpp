/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: collect_data_for_dataset.hpp
 * Project: tmp
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 3rd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Concept of CollectDataForDataset and monoidal structs which can be used by ParallelWriteData
 * structs to collect data within parallel loops and write it to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_TMP_COLLECT_DATA_FOR_DATASET_HPP_
#define LIBS_OBSERVERS2_TMP_COLLECT_DATA_FOR_DATASET_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../kokkosaliases.hpp"
#include "zarr2/dataset.hpp"

/**
 * @brief Concept for CollectDataForDataset is all types that have functions for creating a functor
 * to collect data from a gridbox and/or superdroplets (to use in a Kokkos parallel loop) and then
 * write the data to arrays in a dataset.
 *
 * @tparam CDD The type that satisfies the CollectDataForDataset concept.
 */
template <typename CDD, typename Store>
concept CollectDataForDataset =
    requires(CDD cdd, const Dataset<Store> &ds, const viewd_constgbx d_gbxs,
             const viewd_constsupers totsupers) {
      { cdd.get_functor(d_gbxs, totsupers) };
      { cdd.write_to_arrays(ds) } -> std::same_as<void>;
      { cdd.write_arrayshapes(ds) } -> std::same_as<void>;
    };

/**
 * @brief struct is a new CollectDataForDataset formed from the combination of two structs that
 * also satisfy the CollectDataForDataset concept given the same Store type. Struct that does the
 * actions of the original structs in sequence.
 *
 * @tparam CollectData1 The type of the first CollectDataForDataset.
 * @tparam CollectData2 The type of the second CollectDataForDataset.
 */
template <typename Store, CollectDataForDataset<Store> CollectData1,
          CollectDataForDataset<Store> CollectData2>
struct CombinedCollectDataForDataset {  // TODO(CB) generalise
 private:
  CollectData1 a; /**< The first instance of type of CollectDataForDataset. */
  CollectData2 b; /**< The second instance of type of CollectDataForDataset. */

 public:
  struct Functor {
    CollectData1::Functor a_functor;
    CollectData2::Functor b_functor;

    explicit Functor(const CollectData1 a, const CollectData2 b, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers)
        : a_functor(a.get_functor(d_gbxs, totsupers)),
          b_functor(b.get_functor(d_gbxs, totsupers)) {}

    /* Functor operator to perform copy of each element in parallel in Kokkos Range Policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t nn) const {
      a_functor(nn);
      b_functor(nn);
    }

    /* Functor operator to perform copy of each element in parallel in Kokkos Team Policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {
      a_functor(team_member);
      b_functor(team_member);
    }
  };

  /**
   * @brief Constructs a CombinedCollectDataForDataset object.
   *
   * @param a The first type of CollectDataForDataset instance.
   * @param b The second type of CollectDataForDataset instance.
   */
  CombinedCollectDataForDataset(const CollectData1 a, const CollectData2 b) : a(a), b(b) {}

  Functor get_functor(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    return Functor(a, b, d_gbxs, totsupers);
  }

  void write_to_arrays(const Dataset<Store> &dataset) const {
    a.write_to_arrays(dataset);
    b.write_to_arrays(dataset);
  }

  void write_arrayshapes(const Dataset<Store> &dataset) const {
    a.write_arrayshapes(dataset);
    b.write_arrayshapes(dataset);
  }
};

#endif  // LIBS_OBSERVERS2_TMP_COLLECT_DATA_FOR_DATASET_HPP_
