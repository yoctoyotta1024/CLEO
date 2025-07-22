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
#include "observers/collect_data_for_dataset.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/observers.hpp"
#include "observers/parallel_write_data.hpp"
#include "superdrops/sdmmonitor.hpp"

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
  void before_timestepping(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
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
   * @param d_supers View of superdroplets.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const subviewd_constsupers d_supers) const {
    parallel_write(d_gbxs, d_supers);
  }

  /**
   * @brief Get null monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer that does nothing
   */
  SDMMonitor auto get_sdmmonitor() const { return NullSDMMonitor{}; }
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
  return ConstTstepObserver(interval, DoWriteToDataset(parallel_write));
}

/**
 * @brief Constructs an observer to write data from gridboxes to arrays in a dataset at a
 * constant time interval using a range policy parallelism over the gridboxes.
 * @tparam CollectData Type of collect data function.
 * @param interval Constant timestep interval.
 * @param dataset Dataset to write data to.
 * @param collect_data Function to collect data from gridboxes.
 * @return Constructed observer.
 */
template <typename Dataset, CollectDataForDataset<Dataset> CollectData>
inline Observer auto WriteToDatasetObserver(const unsigned int interval, const Dataset &dataset,
                                            CollectData collect_data) {
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesRangePolicyFunc{}, dataset, collect_data);
  return ConstTstepObserver(interval, DoWriteToDataset(parallel_write));
}

/**
 * @brief Constructs an observer to write data from superdroplets to ragged arrays in a dataset
 * at a constant time interval.
 * @tparam CollectData Type of collect data function.
 * @tparam RaggedCount Type of collect ragged count function.
 * @param interval Constant timestep interval.
 * @param dataset Dataset to write data to.
 * @param collect_data Function to collect data from superdroplets.
 * @param ragged_count Function to collect ragged count for superdroplet data arrays(s).
 * @return Constructed observer.
 */
template <typename Dataset, CollectDataForDataset<Dataset> CollectData,
          CollectRaggedCount<Dataset> RaggedCount>
inline Observer auto WriteToDatasetObserver(const unsigned int interval, const Dataset &dataset,
                                            CollectData collect_data, RaggedCount ragged_count) {
  const auto parallel_write = ParallelWriteSupers(dataset, collect_data, ragged_count);
  return ConstTstepObserver(interval, DoWriteToDataset(parallel_write));
}

#endif  // LIBS_OBSERVERS_WRITE_TO_DATASET_OBSERVER_HPP_
