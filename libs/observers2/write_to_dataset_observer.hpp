/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: write_to_dataset_observer.hpp
 * Project: observers2
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
 * Template for an observer that writes array(s) in a dataset at the start of each step at a
 * constant time interval
 */

#ifndef LIBS_OBSERVERS2_WRITE_TO_DATASET_OBSERVER_HPP_
#define LIBS_OBSERVERS2_WRITE_TO_DATASET_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>

#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "zarr2/dataset.hpp"

/* template class for to call a function during the at_start_step function of an observer
in order to collect variables from gridboxes and/or superdroplets in parallel and then write them
to arrays in a dataset at a constant time interval. */
template <typename Store, typename ParallelWriteData>
class DoWriteInDataset {
 private:
  const Dataset<Store> &dataset;     ///< dataset to write data to
  ParallelWriteData parallel_write;  ///< function like object to call during at_start_step

 public:
  DoWriteInDataset(const Dataset<Store> &dataset, ParallelWriteData parallel_write)
      : dataset(dataset), parallel_write(parallel_write) {}

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes write in dataset observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    parallel_write(dataset, d_gbxs, totsupers);
  }
};

/* constructs observer which writes number of superdrops in each gridbox with a constant
timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store, typename ParallelWriteData>
inline Observer auto WriteInDatasetObserver(const unsigned int interval,
                                            const Dataset<Store> &dataset,
                                            ParallelWriteData parallel_write) {
  const auto obsfunc = DoWriteInDataset(dataset, parallel_write);

  return ConstTstepObserver(interval, obsfunc);
}

#endif  // LIBS_OBSERVERS2_WRITE_TO_DATASET_OBSERVER_HPP_
