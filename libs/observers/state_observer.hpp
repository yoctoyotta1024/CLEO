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
 * Last Modified: Tuesday 9th April 2024
 * Modified By: CB
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
#include "./collect_data_for_dataset.hpp"
#include "./generic_collect_data.hpp"
#include "./thermo_observer.hpp"
#include "./windvel_observer.hpp"
#include "./write_to_dataset_observer.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/dataset.hpp"

/* constructs observer which writes writes thermodynamic and wind velocity variables from the state
of each gridbox with a constant timestep 'interval' using an instance of the WriteToDatasetObserver
class */
template <typename Store>
inline Observer auto StateObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                   const int maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto thermo = CollectThermo(dataset, maxchunk, ngbxs);
  const CollectDataForDataset<Store> auto windvel = CollectWindVel(dataset, maxchunk, ngbxs);

  const CollectDataForDataset<Store> auto collect_data = windvel >> thermo;

  return WriteToDatasetObserver(interval, dataset, collect_data);
}

#endif  // LIBS_OBSERVERS_STATE_OBSERVER_HPP_
