/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: yac_comms.cpp
 * Project: coupldyn_yac
 * Created Date: Friday 3rd May 2024
 * Author: Wilton Loch (WL)
 * Additional Contributors:  Clara Bayley (CB)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * send and receive dynamics functions
 * for SDM when coupled to the yac
 * dynamics solver
 */

#include "coupldyn_yac/yac_comms.hpp"

/* update Gridboxes' states using information received
from YacDynamics solver for 1-way coupling to CLEO SDM.
Kokkos::parallel_for([...]) (on host) is equivalent to:
for (size_t ii(0); ii < ngbxs; ++ii){[...]}
when in serial
// TODO(ALL): make ii indexing compatible with MPI domain decomposition
*/
template <typename GbxMaps, typename CD>
void YacComms::receive_dynamics(const GbxMaps &gbxmaps, const YacDynamics &ffdyn,
                                const viewh_gbx h_gbxs) const {
  const size_t ngbxs(h_gbxs.extent(0));

  Kokkos::parallel_for(
      "receive_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs),
      [=, *this](const size_t ii) { update_gridbox_state(ffdyn, ii, h_gbxs(ii)); });
}

/* updates the state of a gridbox using information
received from YacDynamics solver for 1-way
coupling to CLEO SDM */
void YacComms::update_gridbox_state(const YacDynamics &ffdyn, const size_t ii, Gridbox &gbx) const {
  State &state(gbx.state);

  state.press = ffdyn.get_press(ii);
  state.temp = ffdyn.get_temp(ii);
  state.qvap = ffdyn.get_qvap(ii);
  state.qcond = ffdyn.get_qcond(ii);

  state.wvel = ffdyn.get_wvel(ii);
  state.uvel = ffdyn.get_uvel(ii);
  state.vvel = ffdyn.get_vvel(ii);
}

template void YacComms::send_dynamics<CartesianMaps, YacDynamics>(const CartesianMaps &,
                                                                  const viewh_constgbx,
                                                                  YacDynamics &) const;

template void YacComms::receive_dynamics<CartesianMaps, YacDynamics>(const CartesianMaps &,
                                                                     const YacDynamics &,
                                                                     const viewh_gbx) const;
