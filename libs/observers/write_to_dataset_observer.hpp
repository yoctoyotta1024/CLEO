/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: write_to_dataset_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 8th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Template for an observer that writes array(s) in a dataset at the start of each step at a
 * constant time interval
 */

#ifndef LIBS_OBSERVERS_WRITE_TO_DATASET_OBSERVER_HPP_
#define LIBS_OBSERVERS_WRITE_TO_DATASET_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>

#include "../kokkosaliases.hpp"
#include "./collect_data_for_dataset.hpp"
#include "./const_step_observer.hpp"
#include "./observers.hpp"
#include "./parallel_write_data.hpp"
#include "zarr/dataset.hpp"

/**
 * @class DoWriteToDataset
 * @brief Templated class for writing data from gridboxes and/or superdroplets to a dataset at the
 * constant time intervals by calling the operator of the ParallelWriteData type at the start of
 * each step.
 * @tparam ParallelWriteData Type of function-like object to call during at_start_step.
 */
template <typename ParallelWriteData>
class DoWriteToDataset {
 private:
  ParallelWriteData parallel_write; /**< Function-like object to call during at_start_step. */

 public:
  /**
   * @brief Constructor for DoWriteToDataset.
   * @param parallel_write Function-like object to call during at_start_step.
   */
  explicit DoWriteToDataset(ParallelWriteData parallel_write) : parallel_write(parallel_write) {}

  /**
   * @brief Placeholder for before timestepping functionality and to make class satisfy observer
   * concept.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write in dataset observer\n";
  }

  /**
   * @brief Placeholder for after timestepping functionality and to make class satisfy observer
   * concept.
   */
  void after_timestepping() const {}

  /**
   * @brief Calls the parallel_write function during at_start_step.
   * @param t_mdl Current model time.
   * @param d_gbxs View of gridboxes.
   * @param totsupers View of superdroplets.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    parallel_write(d_gbxs, totsupers);
  }

  /**
   * @brief No operation at the start of a SDM substep.
   *
   * @param t_sdm The unsigned int parameter representing the current model time.
   * @param d_gbxs The view of gridboxes in device memory.
   */
  void at_start_sdm_substep(const unsigned int t_sdm, const viewd_constgbx d_gbxs) const {}
};

/**
 * @brief Constructs an observer to write data from gridboxes and/or superdroplets to a dataset at a
 * constant time interval according to the ParallelWriteData struct.
 * @tparam ParallelWriteData Type of function-like object to call during at_start_step.
 * @param interval Constant timestep interval.
 * @param parallel_write Function-like object to call during at_start_step.
 * @return Constructed observer.
 */
template <typename ParallelWriteData>
inline Observer auto WriteToDatasetObserver(const unsigned int interval,
                                            ParallelWriteData parallel_write) {
  return ConstStepObserver(interval, DoWriteToDataset(parallel_write));
}

/**
 * @brief Constructs an observer to write data from gridboxes to arrays in a dataset at a
 * constant time interval using a range policy parallelism over the gridboxes.
 * @tparam Store Type of store for dataset.
 * @tparam CollectData Type of collect data function.
 * @param interval Constant timestep interval.
 * @param dataset Dataset to write data to.
 * @param collect_data Function to collect data from gridboxes.
 * @return Constructed observer.
 */
template <typename Store, CollectDataForDataset<Store> CollectData>
inline Observer auto WriteToDatasetObserver(const unsigned int interval,
                                            const Dataset<Store> &dataset,
                                            CollectData collect_data) {
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesRangePolicyFunc{}, dataset, collect_data);
  return ConstStepObserver(interval, DoWriteToDataset(parallel_write));
}

/**
 * @brief Constructs an observer to write data from superdroplets to ragged arrays in a dataset
 * at a constant time interval.
 * @tparam Store Type of store for dataset.
 * @tparam CollectData Type of collect data function.
 * @tparam RaggedCount Type of collect ragged count function.
 * @param interval Constant timestep interval.
 * @param dataset Dataset to write data to.
 * @param collect_data Function to collect data from superdroplets.
 * @param ragged_count Function to collect ragged count for superdroplet data arrays(s).
 * @return Constructed observer.
 */
template <typename Store, CollectDataForDataset<Store> CollectData,
          CollectRaggedCount<Store> RaggedCount>
inline Observer auto WriteToDatasetObserver(const unsigned int interval,
                                            const Dataset<Store> &dataset, CollectData collect_data,
                                            RaggedCount ragged_count) {
  const auto parallel_write = ParallelWriteSupers(dataset, collect_data, ragged_count);
  return ConstStepObserver(interval, DoWriteToDataset(parallel_write));
}

#endif  // LIBS_OBSERVERS_WRITE_TO_DATASET_OBSERVER_HPP_
