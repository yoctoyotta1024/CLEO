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

#ifndef LIBS_OBSERVERS2_STATE_OBSERVER_HPP_
#define LIBS_OBSERVERS2_STATE_OBSERVER_HPP_

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
#include "./write_to_dataset_observer.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/* Operator is functor to perform copy of pressure in each gridbox to d_data in parallel.
Note conversion of pressure from double (8 bytes) to single precision (4 bytes float) in output */
struct PressFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto press = static_cast<float>(d_gbxs(ii).state.press);
    d_data(ii) = press;
  }
};

/* returns CollectDataForDataset which writes a state variable from
each gridbox to an array in a dataset in a given store for a given datatype and using a given
function-like functor */
template <typename Store, typename T, typename FunctorFunc>
CollectDataForDataset<Store> auto CollectStateVariable(
    const Dataset<Store> &dataset, const FunctorFunc ffunc, const std::string_view name,
    const std::string_view units, const std::string_view dtype, const double scale_factor,
    const size_t maxchunk, const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr =
      dataset.template create_array<T>(name, units, dtype, scale_factor, chunkshape, dimnames);
  return GenericCollectData(ffunc, xzarr, ngbxs);
}

/* constructs observer which writes writes thermodynamic variables from the state of each gridbox
with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
template <typename Store>
inline Observer auto ThermoObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                    const int maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto press = CollectStateVariable<Store, float, PressFunc>(
      dataset, PressFunc{}, "press", "hPa", "<f4", dlc::P0 / 100, maxchunk, ngbxs);

  const auto collect_thermodata = press;  // TODO(CB) WIP
  return WriteToDatasetObserver(interval, dataset, collect_thermodata);
}

// /* constructs observer which writes writes the wind velocity from the state of each gridbox
// with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
// template <typename Store>
// inline Observer auto WindObserver(const unsigned int interval, const Dataset<Store> &dataset,
//                                   const int maxchunk, const size_t ngbxs) {
//   // TODO(CB) WIP
//   return WriteToDatasetObserver(interval, dataset, collect_winddata);
// }

// /* constructs observer which writes writes thermodynamic and wind velocity variables from the
// state of each gridbox with a constant timestep 'interval' using an instance of the
// WriteToDatasetObserver class */ template <typename Store> inline Observer auto
// StateObserver(const unsigned int interval, const Dataset<Store> &dataset,
//                                    const int maxchunk, const size_t ngbxs) {
//   // TODO(CB) WIP
//   return WriteToDatasetObserver(interval, dataset, collect_statedata);
// }

#endif  // LIBS_OBSERVERS2_STATE_OBSERVER_HPP_
