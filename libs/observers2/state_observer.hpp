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
 * Last Modified: Thursday 28th March 2024
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

#include <Kokkos_Core.hpp>
#include <cassert>
#include <concepts>
#include <iostream>
#include <memory>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "./write_gridboxes.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"

// returns GridboxDataWriter which writes the pressure in each gridbox to an array in a dataset in a
// store
template <typename Store>
GridboxDataWriter<Store> auto PressureGbxWriter(Dataset<Store> &dataset, const int maxchunk,
                                                const size_t ngbxs) {
  // Operator is functor to perform copy of pressure in each gridbox to d_data in parallel
  struct PressureFunc {
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii, viewd_constgbx d_gbxs,
                    Buffer<double>::mirrorviewd_buffer d_data) const {
      d_data(ii) = d_gbxs(ii).state.press;
    }
  };

  // create shared pointer to 2-D array in dataset for pressure in each gridbox over time
  std::shared_ptr<XarrayZarrArray<Store, double>> xzarr_ptr =
      std::make_shared<XarrayZarrArray<Store, double>>(dataset.template create_array<double>(
          "press", "hPa", "<f8", dlc::P0 / 100, good2Dchunkshape(maxchunk, ngbxs),
          {"time", "gbxindex"}));

  return GenericGbxWriter<Store, double, PressureFunc>(dataset, PressureFunc{}, xzarr_ptr, ngbxs);
}

/* constructs observer which writes variables from the state of each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto StateObserver(const unsigned int interval, Dataset<Store> &dataset,
                                   const int maxchunk, const size_t ngbxs) {
  const auto c = CombineGDW<Store>{};

  auto pressure_writer = get_pressure_gbxwriter(dataset, maxchunk, ngbxs);
  const GridboxDataWriter<Store> auto writer = c(pressure_writer, NullGbxWriter<Store>{});

  return ConstTstepObserver(interval, WriteGridboxes(dataset, writer));
}

#endif  // LIBS_OBSERVERS2_STATE_OBSERVER_HPP_
