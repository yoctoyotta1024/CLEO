/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: yaccomms.hpp
 * Project: coupldyn_yac
 * Created Date: Tuesday 31st October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * send and receive dynamics functions
 * for SDM when coupled to the yac
 * dynamics solver
 */

#ifndef LIBS_COUPLDYN_YAC_YAC_COMMS_HPP_
#define LIBS_COUPLDYN_YAC_YAC_COMMS_HPP_

#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "coupldyn_yac/yac_cartesian_dynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"

/* 1-way coupling from coupldyn to CLEO's gridboxes where
coupldyn sends information to h_gbxs but doesn't
receive any back. Struct obeys coupling comms concept */
struct YacComms {
 private:
  /* updates the state of a gridbox using information
  received from YacCartesianDynamics solver for 1-way
  coupling to CLEO SDM */
  void update_gridbox_state(const YacCartesianDynamics &ffdyn, const size_t ii, Gridbox &gbx) const;

 public:
  /* send information from Gridboxes' states
  to coupldyn is null for YacCartesianDynamics*/
  template <typename GbxMaps, typename CD = YacCartesianDynamics>
  void send_dynamics(const GbxMaps &, const viewh_constgbx,
                     const YacCartesianDynamics &ffdyn) const;

  /* update Gridboxes' states using information
  received from YacCartesianDynamics solver for
  1-way coupling to CLEO SDM */
  template <typename GbxMaps, typename CD = YacCartesianDynamics>
  void receive_dynamics(const GbxMaps &gbxmaps, const YacCartesianDynamics &ffdyn,
                        const viewh_gbx h_gbxs) const;
};

#endif  // LIBS_COUPLDYN_YAC_YAC_COMMS_HPP_
