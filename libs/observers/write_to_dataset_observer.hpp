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
 * Last Modified: Tuesday 9th April 2024
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
#include "./observers.hpp"
#include "./parallel_write_data.hpp"
#include "zarr/dataset.hpp"

/* template class for to call a function during the at_start_step function of an observer
in order to collect variables from gridboxes and/or superdroplets in parallel and then write them
to arrays in a dataset at a constant time interval. */
template <typename ParallelWriteData>
class DoWriteToDataset {
 private:
  ParallelWriteData parallel_write;  ///< function like object to call during at_start_step

 public:
  explicit DoWriteToDataset(ParallelWriteData parallel_write) : parallel_write(parallel_write) {}

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write in dataset observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    parallel_write(d_gbxs, totsupers);
  }
};

/* constructs observer which writes some data in parallel according to parallel_write at the start
of a constant timestep 'interval' using instances of the ConstTstepObserver and DoWriteToDataset
classes */
template <typename ParallelWriteData>
inline Observer auto WriteToDatasetObserver(const unsigned int interval,
                                            ParallelWriteData parallel_write) {
  return ConstTstepObserver(interval, DoWriteToDataset(parallel_write));
}

/* constructs observer which writes some data from gridboxes in parallel into arrays in a dataset
according to parallel_write at the start of a constant timestep 'interval' using instances of the
ConstTstepObserver and DoWriteToDataset classes */
template <typename Store, CollectDataForDataset<Store> CollectData>
inline Observer auto WriteToDatasetObserver(const unsigned int interval,
                                            const Dataset<Store> &dataset,
                                            CollectData collect_data) {
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesRangePolicyFunc{}, dataset, collect_data);
  return ConstTstepObserver(interval, DoWriteToDataset(parallel_write));
}

/* constructs observer which writes some data from superdroplets in parallel into ragged arrays in a
dataset according to parallel_write at the start of a constant timestep 'interval' using instances
of the ConstTstepObserver and DoWriteToDataset classes */
template <typename Store, CollectDataForDataset<Store> CollectData,
          CollectRaggedCount<Store> RaggedCount>
inline Observer auto WriteToDatasetObserver(const unsigned int interval,
                                            const Dataset<Store> &dataset, CollectData collect_data,
                                            RaggedCount ragged_count) {
  const auto parallel_write = ParallelWriteSupers(dataset, collect_data, ragged_count);
  return ConstTstepObserver(interval, DoWriteToDataset(parallel_write));
}

#endif  // LIBS_OBSERVERS_WRITE_TO_DATASET_OBSERVER_HPP_
