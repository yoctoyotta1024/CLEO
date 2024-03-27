/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state_obs.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 27th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables related to Gridboxes' state at the start of
 * each timestep to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_STATE_OBS_HPP_
#define LIBS_OBSERVERS2_STATE_OBS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"
#include "zarr2/dataset.hpp"

Buffer<double>::viewh_buffer copy_press(const viewd_constgbx d_gbxs) {
  const auto ngbxs = size_t{d_gbxs.extent(0)};
  auto h_press = Buffer<double>::viewh_buffer("press", ngbxs);
  auto d_press = Kokkos::create_mirror_view(ExecSpace(), h_press);

  Kokkos::parallel_for(
      "stateobs", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
      KOKKOS_LAMBDA(const size_t ii) { d_press(ii) = d_gbxs(ii).state.press; });
  Kokkos::deep_copy(h_press, d_press);

  return h_press;
}

/* observe variables in the state of each
gridbox and write them to repspective arrays
in a store as determined by the Dataset */
template <typename Store>
class DoStateObs {
 private:
  Dataset<Store> &dataset;
  XarrayZarrArray<Store, double> xzarr_press;

 public:
  DoStateObs(Dataset<Store> &dataset, const std::vector<size_t> &chunkshape)
      : dataset(dataset),
        xzarr_press(dataset.template create_array<double>("press", "hPa", "<f8", dlc::P0 / 100,
                                                          chunkshape, {"time", "gbxindex"})) {}

  ~DoStateObs() { dataset.write_arrayshape(xzarr_press); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes State observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    // TODO(CB) complete function WIP
    const auto h_press = copy_press(d_gbxs);
    dataset.write_to_array(xzarr_press, h_press);

    // dataset.set_dimension({"time", time+1});
    // dataset.write_arrayshape(xzarr_press);
  }
};

/* constructs observer of variables in the state
of each gridbox with a constant timestep 'interval'
using an instance of the DoStateObs class */
template <typename Store>
inline Observer auto StateObserver(const unsigned int interval, Dataset<Store> &dataset,
                                   const int maxchunk, const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto obs = DoStateObs<Store>(dataset, chunkshape);
  return ConstTstepObserver(interval, obs);
}

#endif  // LIBS_OBSERVERS2_STATE_OBS_HPP_
