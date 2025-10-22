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
  const int coupling_flag = ffdyn.get_dynvars()->get_yac_coupling_flag();

  if ((coupling_flag ==2) || (coupling_flag ==1)) {
  ffdyn.get_dynvars()->receive_fields_from_yac();
  }
  Kokkos::parallel_for(
      "receive_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs),
      [=, *this](const size_t ii) { update_gridbox_state(ffdyn, ii, h_gbxs(ii)); });
}

template <typename GbxMaps, typename CD>
void YacComms::send_dynamics(const GbxMaps &gbxmaps, const viewh_constgbx h_gbxs,
    const YacDynamics &ffdyn) const {
  const size_t ngbxs(h_gbxs.extent(0));

  Kokkos::View<double*, HostSpace> temp_state("temp_state", ngbxs);
  Kokkos::View<double*, HostSpace> qvap_state("qvap_state", ngbxs);
  Kokkos::View<double*, HostSpace> qcond_state("qcond_state", ngbxs);

  Kokkos::parallel_for(
      "send_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs),
      [=, *this](const size_t ii) {
      State &state(h_gbxs(ii).state);
      temp_state[ii] = state.temp;
      qvap_state[ii] = state.qvap;
      qcond_state[ii] = state.qcond;
      });
  const int coupling_flag = ffdyn.get_dynvars()->get_yac_coupling_flag();

  if (coupling_flag ==2) {
  ffdyn.get_dynvars()->send_fields_to_yac(temp_state.data(), qvap_state.data(), qcond_state.data());
  }
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
                                                                  const YacDynamics &) const;

template void YacComms::receive_dynamics<CartesianMaps, YacDynamics>(const CartesianMaps &,
                                                                     const YacDynamics &,
                                                                     const viewh_gbx) const;
