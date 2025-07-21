/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to write variables related to Gridboxes' state at the start of
 * a constant interval timestep to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS_STATE_OBSERVER_HPP_
#define LIBS_OBSERVERS_STATE_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "observers/collect_data_for_dataset.hpp"
#include "observers/generic_collect_data.hpp"
#include "observers/thermo_observer.hpp"
#include "observers/windvel_observer.hpp"
#include "observers/write_to_dataset_observer.hpp"

/**
 * @brief Constructs an observer which writes the state of a gridbox (thermodynamics and
 * wind velocity components) in each gridbox at start of each observation timestep to an array with
 * a constant observation timestep "interval".
 *
 * This function collects thermodynamic properties and wind velocities from the dataset and combines
 * them into a single collection of state data.
 *
 * @tparam Dataset Type of dataset
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Observer An observer instance for writing the state data.
 */
template <typename Dataset>
inline Observer auto StateObserver(const unsigned int interval, const Dataset &dataset,
                                   const size_t maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Dataset> auto thermo = CollectThermo(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Dataset> auto windvel = CollectWindVel(dataset, maxchunk, ngbxs);

  const CollectDataForDataset<Dataset> auto collect_data =
      CombinedCollectDataForDataset(windvel, thermo);

  return WriteToDatasetObserver(interval, dataset, collect_data);
}

#endif  // LIBS_OBSERVERS_STATE_OBSERVER_HPP_
