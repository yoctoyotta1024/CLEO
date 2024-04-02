/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state_observer.hpp
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
 * Observer to output variables related to Gridboxes' state at the start of
 * each timestep to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_STATE_OBSERVER_HPP_
#define LIBS_OBSERVERS2_STATE_OBSERVER_HPP_

#include <concepts>

#include "./observers.hpp"
#include "./state_writers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "./write_gridboxes.hpp"
#include "zarr2/dataset.hpp"

/* constructs observer which writes thermodynamic variables from the state of each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto ThermoObserver(const unsigned int interval, Dataset<Store> &dataset,
                                    const int maxchunk, const size_t ngbxs) {
  const WriteGridboxToArray<Store> auto thermowriter = ThermoWriter(dataset, maxchunk, ngbxs);
  const auto looppolicy = ParallelGbxsRangePolicy{};
  return ConstTstepObserver(interval, WriteGridboxes(dataset, thermowriter, looppolicy));
}

/* constructs observer which writes the wind velocity from the state of each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto WindObserver(const unsigned int interval, Dataset<Store> &dataset,
                                  const int maxchunk, const size_t ngbxs) {
  const WriteGridboxToArray<Store> auto windwriter = WindVelocityWriter(dataset, maxchunk, ngbxs);
  const auto looppolicy = ParallelGbxsRangePolicy{};
  return ConstTstepObserver(interval, WriteGridboxes(dataset, windwriter, looppolicy));
}

/* constructs observer which writes variables from the state of each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto StateObserver(const unsigned int interval, Dataset<Store> &dataset,
                                   const int maxchunk, const size_t ngbxs) {
  const WriteGridboxToArray<Store> auto thermowriter = ThermoWriter(dataset, maxchunk, ngbxs);
  const WriteGridboxToArray<Store> auto windwriter = WindVelocityWriter(dataset, maxchunk, ngbxs);
  const WriteGridboxToArray<Store> auto statewriter = CombineGDW<Store>{}(thermowriter, windwriter);
  const auto looppolicy = ParallelGbxsRangePolicy{};

  return ConstTstepObserver(interval, WriteGridboxes(dataset, statewriter, looppolicy));
}

#endif  // LIBS_OBSERVERS2_STATE_OBSERVER_HPP_
