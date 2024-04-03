/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state_observer.hpp
 * Project: tmp
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 4th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to write variables related to Gridboxes' state at the start of
 * a constant interval timestep to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_TMP_STATE_OBSERVER_HPP_
#define LIBS_OBSERVERS2_TMP_STATE_OBSERVER_HPP_

#include "./write_to_dataset_observer.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/* constructs observer which writes writes thermodynamic variables from the state of each gridbox
with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
template <typename Store>
inline Observer auto ThermoObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                    const int maxchunk, const size_t ngbxs) {
  return WriteToDatasetObserver(interval, dataset, collect_thermodata);
}

/* constructs observer which writes writes the wind velocity from the state of each gridbox
with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
template <typename Store>
inline Observer auto WindObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                  const int maxchunk, const size_t ngbxs) {
  return WriteToDatasetObserver(interval, dataset, collect_winddata);
}

/* constructs observer which writes writes thermodynamic and wind velocity variables from the state
of each gridbox with a constant timestep 'interval' using an instance of the WriteToDatasetObserver
class */
template <typename Store>
inline Observer auto StateObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                   const int maxchunk, const size_t ngbxs) {
  return WriteToDatasetObserver(interval, dataset, collect_statedata);
}

#endif  // LIBS_OBSERVERS2_TMP_STATE_OBSERVER_HPP_
