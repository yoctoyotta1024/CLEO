/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: massmoments_observer.hpp
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
 * Observer to output the mass moments of the droplet size distribution in each gridbox
 * to individual arrays in a dataset a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS2_TMP_MASSMOMENTS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_TMP_MASSMOMENTS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "./observers.hpp"
#include "./write_to_dataset_observer.hpp"
#include "zarr2/dataset.hpp"

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
template <typename Store>
inline Observer auto MassMomentsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                         const int maxchunk, const size_t ngbxs) {
  return WriteInDatasetObserver(interval, dataset, write_massmoments);
}

/* constructs observer which writes mass moments of raindroplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the WritetoDatasetObserver class */
template <typename Store>
inline Observer auto MassMomentsRaindropsObserver(const unsigned int interval,
                                                  const Dataset<Store> &dataset, const int maxchunk,
                                                  const size_t ngbxs) {
  return WriteInDatasetObserver(interval, dataset, write_massmoments_raindrops);
}

#endif  // LIBS_OBSERVERS2_TMP_MASSMOMENTS_OBSERVER_HPP_
