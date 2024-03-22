/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: stateobs.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 22nd March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observers to output variables from gridboxes' state to arrays in a zarr file system storage
 */

#ifndef LIBS_OBSERVERS_STATEOBS_HPP_
#define LIBS_OBSERVERS_STATEOBS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"
#include "zarr/statebuffers.hpp"
#include "zarr/twodstorage.hpp"
#include "zarr/windbuffers.hpp"

/* constructs observer of variables in the state
of each gridbox with a constant timestep 'interval'
using an instance of the DoStateObs class */
inline Observer auto StateObserver(const unsigned int interval, FSStore &store, const int maxchunk,
                                   const size_t ngbxs);

/* observe variables in the state of each
gridbox and write them to repspective arrays
in a store as determined by the StateBuffers
and TwoDMulitVarStorage types */
class DoStateObs {
 private:
  using store_type_thermo = TwoDMultiVarStorage<StateBuffers<double>, State>;
  using store_type_winds = TwoDMultiVarStorage<WindBuffers<double>, State>;
  std::shared_ptr<store_type_thermo> zarr_thermo;
  std::shared_ptr<store_type_winds> zarr_winds;

 public:
  DoStateObs(FSStore &store, const int maxchunk, const size_t ngbxs)
      : zarr_thermo(std::make_shared<store_type_thermo>(store, maxchunk, "<f8", ngbxs, "")),
      : zarr_winds(std::make_shared<store_type_winds>(store, maxchunk, "<f8", ngbxs, "")) {}

  void before_timestepping(const viewh_constgbx h_gbxs) const {
    std::cout << "observer includes StateObserver\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {}

  /* writes some variables from gridbox state
  to 2-D zarr storages as determined by the
  StateBuffers struct */
  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {
    zarr_thermo->values_to_storage(gbx.state);
    zarr_winds->values_to_storage(gbx.state);
  }
};

/* constructs observer of variables in the state
of each gridbox with a constant timestep 'interval'
using an instance of the DoStateObs class */
inline Observer auto StateObserver(const unsigned int interval, FSStore &store, const int maxchunk,
                                   const size_t ngbxs) {
  const auto obs = DoStateObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

#endif  // LIBS_OBSERVERS_STATEOBS_HPP_
