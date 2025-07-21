/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: parallel_write_data.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Some "ParallelWriteData" function-like objects (see write_to_dataset_observer.hpp) for
 * writing data from gridboxes and/or superdroplets to arrays in a dataset.
 */

#ifndef LIBS_OBSERVERS_PARALLEL_WRITE_DATA_HPP_
#define LIBS_OBSERVERS_PARALLEL_WRITE_DATA_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../kokkosaliases.hpp"

namespace KCS = KokkosCleoSettings;

/**
 * @brief Struct for a function-like object with operator() to call when using type as
 * parallel_gridboxes_func object in ParallelWriteGridboxes.
 *
 * Struct for a function-like object with operator() suitable for parallel_gridboxes_func object in
 * ParallelWriteGridboxes in order to loop over gridbxoes using a kokkos range policy.
 */
struct ParallelGridboxesRangePolicyFunc {
  /**
   * @brief Parallel loop over gridboxes using a Kokkos Range Policy.
   *
   * Kokkos::parallel_for([...]) is equivalent in serial to:
   * for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}.
   *
   *  _Note:_ type for Functor used in this call must have operator() with signature:
   * void operator()(const size_t ii).
   *
   * @tparam Functor The type of the functor.
   * @param functor The functor to be executed in parallel.
   * @param d_gbxs The view of gridboxes.
   */
  template <typename Functor>
  void operator()(const Functor functor, const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for("write_gridboxes_range", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
                         functor);
  }
};

/**
 * @brief Struct for a function-like object with operator() to call when using type as
 * parallel_gridboxes_func object in ParallelWriteGridboxes.
 *
 * Struct for a function-like object with operator() suitable for parallel_gridboxes_func object in
 * ParallelWriteGridboxes in order to loop over gridbxoes using a kokkos team policy.
 */
struct ParallelGridboxesTeamPolicyFunc {
  /**
   * @brief Parallel loop over gridboxes using a Kokkos Team Policy.
   *
   * Kokkos::parallel_for([...]) is equivalent in serial to:
   * for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}.
   *
   *  _Note:_ type for Functor used in this call must have operator() with signature:
   * void operator()(const TeamMember &team_member).
   *
   * @tparam Functor The type of the functor.
   * @param functor The functor to be executed in parallel.
   * @param d_gbxs The view of gridboxes.
   */
  template <typename Functor>
  void operator()(const Functor functor, const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for("write_gridboxes_team", TeamPolicy(ngbxs, KCS::team_size), functor);
  }
};

/**
 * @brief Struct for "ParallelWriteData" (see write_to_dataset_observer.hpp) to collect data from
 * gridboxes in a loop (e.g. using Kokkos::parallel_for with a range or team policy) and then write
 * that data to arrays in a dataset.
 *
 * This struct is responsible for collecting data from gridboxes and writing it to arrays in a
 * dataset.
 *
 * The ParallelGridboxesFunc type is a function-like object responsible for looping over
 * gridboxes in parallel (see ParallelGridboxesRangePolicyFunc or ParallelGridboxesTeamPolicyFunc).
 * The CollectData type satisfies the concept for CollectDataForDataset and the signature for the
 * operator of the type it returns from its get_functor() call should be compatible with the
 * signature of the functor required by the ParallelGridboxesFunc type.
 *
 * @tparam Dataset The type dataset used to write data to a store.
 * @tparam ParallelGridboxesFunc Function-like object for call to loop over gridboxes.
 * @tparam CollectData Object satisfying the CollectDataForDataset concept for the given
 * dataset.
 */
template <typename Dataset, typename ParallelGridboxesFunc,
          CollectDataForDataset<Dataset> CollectData>
class ParallelWriteGridboxes {
 private:
  ParallelGridboxesFunc
      parallel_gridboxes_func; /**< Function-like object for call to loop over gridboxes.*/
  const Dataset &dataset;      /**< Dataset to write data to. */
  CollectData collect_data; /**< CollectData Object satisfying the CollectDataForDataset concept. */

 public:
  /**
   * @brief Constructs a new ParallelWriteGridboxes object.
   *
   * @param parallel_gridboxes_func Function-like object for call to loop over gridboxes.
   * @param dataset The dataset to write data to.
   * @param collect_data The object satisfying the CollectDataForDataset concept to collect data.
   */
  ParallelWriteGridboxes(ParallelGridboxesFunc parallel_gridboxes_func, const Dataset &dataset,
                         CollectData collect_data)
      : parallel_gridboxes_func(parallel_gridboxes_func),
        dataset(dataset),
        collect_data(collect_data) {}

  /**
   * @brief Destructor for the ParallelWriteGridboxes class.
   */
  ~ParallelWriteGridboxes() { collect_data.write_arrayshapes(dataset); }

  /**
   * @brief Executes the operation to collect data from gridboxes and write it to arrays in the
   * dataset.
   *
   * Use the funtor returned by CollectData's get_functor() call to collect data from gridboxes in a
   * parallel loop as determined by the parallel_gridboxes_func operator().
   *
   * Inclusion of d_supers in function signature so that object can be used as "ParallelWriteData"
   * function-like object in DoWriteToDataset struct.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   * @param d_supers The view of superdroplets in device memory.
   */
  void operator()(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    auto functor = collect_data.get_functor(d_gbxs, d_supers);
    parallel_gridboxes_func(functor, d_gbxs);
    collect_data.write_to_arrays(dataset);
  }
};

/**
 * @brief Concept for CollectRaggedCount is all types that have functions for writing the ragged
 * count of superdroplet arrays to an array in a dataset.
 *
 * @tparam CRC The type that satisfies the CollectRaggedCount concept.
 */
template <typename CRC, typename Dataset>
concept CollectRaggedCount =
    requires(CRC crc, const Dataset &ds, const subviewd_constsupers d_supers) {
      { crc.write_to_array(ds, d_supers) } -> std::same_as<void>;
      { crc.write_arrayshape(ds) } -> std::same_as<void>;
    };

/**
 * @brief Struct for "ParallelWriteData" (see write_to_dataset_observer.hpp) to collect data from
 * superdroplets in a (parallel) loop and then write that data to ragged arrays in a dataset.
 *
 * This struct is responsible for collecting data from superdroplets and writing it to ragged arrays
 * in a dataset with a given store.
 *
 * The CollectData type satisfies the concept for CollectDataForDataset and the signature for
 * the operator of the type it returns from its get_functor() call should be compatible with the
 * signature of the functor required by the Kokkos::parallel_for loop over superdroplets.
 *
 * @tparam Dataset The type dataset used to write data to a store.
 * @tparam CollectData The object for collecting data satsifying the CollectDataForDataset concept.
 * @tparam RaggedCount The type of the function object for writing the ragged count variable in the
 * dataset satisfying the CollectRaggedCount concept.
 */
template <typename Dataset, CollectDataForDataset<Dataset> CollectData,
          CollectRaggedCount<Dataset> RaggedCount>
class ParallelWriteSupers {
 private:
  const Dataset &dataset; /**< dataset to write data to */
  CollectData collect_data;
  /**< functions to collect data within loop over superdroplets and write to ragged array(s) */
  RaggedCount ragged_count; /**< functions to write ragged count variable to a dataset */

  /**
   * @brief Parallel loop over superdroplets using a Kokkos Range Policy.
   *
   * Kokkos::parallel_for([...]) is equivalent in serial to:
   * for (size_t kk(0); kk < d_supers.extent(0); ++kk){[...]}
   *
   * _Note:_ type for Functor used in this call must have operator() with signature:
   * void operator()(const size_t kk).
   *
   * @tparam Functor The type of the functor.
   * @param functor The functor to be executed in parallel.
   * @param d_supers The view of superdroplets on device.
   */
  template <typename Functor>
  void parallel_supers_func(const Functor functor, const subviewd_constsupers d_supers) const {
    const size_t totnsupers(d_supers.extent(0));
    Kokkos::parallel_for("write_supers", Kokkos::RangePolicy<ExecSpace>(0, totnsupers), functor);
  }

 public:
  /**
   * @brief Constructs a new ParallelWriteSupers object.
   *
   * @param dataset The dataset to write data to.
   * @param collect_data Object for collecting data satsifying the CollectDataForDataset concept.
   * @param ragged_count Object for writing the ragged count variable in the dataset satisfying the
   * CollectRaggedCount concept.
   */
  ParallelWriteSupers(const Dataset &dataset, CollectData collect_data, RaggedCount ragged_count)
      : dataset(dataset), collect_data(collect_data), ragged_count(ragged_count) {}

  /**
   * @brief Destructor for the ParallelWriteSupers class.
   *
   */
  ~ParallelWriteSupers() {
    collect_data.write_ragged_arrayshapes(dataset);
    ragged_count.write_arrayshape(dataset);
  }

  /* Use the CollectData instance's functor to collect data from superdroplets in a parallel loop.
   */

  /**
   * @brief Executes the operation to collect data from superdroplets and write it to ragged arrays
   * in the dataset.
   *
   * Use the funtor returned by CollectData's get_functor() call to collect data from superdroplets
   * in a parallel loop as determined by the parallel_supers_func function call. Then write the data
   * in the dataset alongside ragged count for arrays.
   *
   * Inclusion of d_gbxs in function signature so that object can be used as "ParallelWriteData"
   * function-like object in DoWriteToDataset struct.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   * @param d_supers The view of superdroplets in device memory.
   */
  void operator()(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    collect_data.reallocate_views(d_supers.extent(0));
    auto functor = collect_data.get_functor(d_gbxs, d_supers);
    parallel_supers_func(functor, d_supers);
    collect_data.write_to_ragged_arrays(dataset);
    ragged_count.write_to_array(dataset, d_supers);
  }
};

#endif  // LIBS_OBSERVERS_PARALLEL_WRITE_DATA_HPP_
