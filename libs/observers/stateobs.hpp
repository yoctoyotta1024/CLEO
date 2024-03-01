/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: stateobs.hpp
 * Project: observers
 * Created Date: Monday 23rd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 23rd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observers to output variables from
 * gridboxes' state to arrays in a zarr
 * file system storage
 */

#ifndef LIBS_OBSERVERS_STATEOBS_HPP_
#define LIBS_OBSERVERS_STATEOBS_HPP_

#include <concepts>
#include <iostream>
#include <memory>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"
#include "zarr/statebuffers.hpp"
#include "zarr/twodstorage.hpp"

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
  using store_type = TwoDMultiVarStorage<StateBuffers<double>, State>;
  std::shared_ptr<store_type> zarr;

 public:
  DoStateObs(FSStore &store, const int maxchunk, const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk, "<f8", ngbxs, "")) {}

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
    zarr->values_to_storage(gbx.state);
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
