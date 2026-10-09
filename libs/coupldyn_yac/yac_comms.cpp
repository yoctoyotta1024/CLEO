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
from YacCartesianDynamics solver for 1-way coupling to CLEO SDM.
Kokkos::parallel_for([...]) (on host) is equivalent to:
for (size_t ii(0); ii < ngbxs; ++ii){[...]}
when in serial
// TODO(ALL): make ii indexing compatible with MPI domain decomposition
*/

template <typename GbxMaps, typename CD>
void YacComms::receive_dynamics(const GbxMaps& gbxmaps, const YacCartesianDynamics& ffdyn,
                                const viewh_gbx h_gbxs) const {
  const int coupling_flag = ffdyn.get_yac_coupling_flag();

  if ((coupling_flag == 2) || (coupling_flag == 1)) {
    ffdyn.receive_fields_from_yac();

    const size_t ngbxs(h_gbxs.extent(0));
    Kokkos::parallel_for(
        "receive_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs),
        [=, *this](const size_t ii) { update_gridbox_state(ffdyn, ii, h_gbxs(ii)); });
  }
}

template <typename GbxMaps, typename CD>
void YacComms::send_dynamics(const GbxMaps& gbxmaps, const viewh_constgbx h_gbxs,
                             const YacCartesianDynamics& ffdyn) const {
  const int coupling_flag = ffdyn.get_yac_coupling_flag();

  if (coupling_flag == 2) {
    const size_t ngbxs(h_gbxs.extent(0));
    Kokkos::parallel_for("send_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs),
                         [=, *this](const size_t ii) {
                           const State& state(h_gbxs(ii).state);
                           const auto temp = state.temp;
                           const auto qvap = state.qvap;
                           const auto qcloud = 0.0;  // total mass of "cloud"
                           const auto qrain = 0.0;   // total mass of "rain"
                           update_send_buffers(ffdyn, ii, temp, qvap, qcloud, qrain);
                         });
    ffdyn.send_fields_to_yac();
  }
}

/* updates the state of a gridbox using information received from YacCartesianDynamics solver
for 1-way or 2-way coupling to from Dynamics to CLEO SDM */
void YacComms::update_gridbox_state(const YacCartesianDynamics& ffdyn, const size_t ii,
                                    Gridbox& gbx) const {
  State& state(gbx.state);

  state.press = ffdyn.get_press(ii);
  state.temp = ffdyn.get_temp(ii);
  state.qvap = ffdyn.get_qvap(ii);

  state.wvel = ffdyn.get_wvel(ii);
  state.uvel = ffdyn.get_uvel(ii);
  state.vvel = ffdyn.get_vvel(ii);
}

/* fills YacComms send buffers with the state of a gridbox for a 2-way
coupling to CLEO SDM */
void YacComms::update_send_buffers(const YacCartesianDynamics& ffdyn, const size_t ii,
                                   const double temp, const double qvap, const double qcloud,
                                   const double qrain) const {
  ffdyn.set_temp_delta(ii, temp - ffdyn.get_temp(ii));
  ffdyn.set_qvap_delta(ii, qvap - ffdyn.get_qvap(ii));
  ffdyn.set_qcloud_delta(ii, qcloud - ffdyn.get_qcloud(ii));
  ffdyn.set_qrain_delta(ii, qrain - ffdyn.get_qrain(ii));
}

template void YacComms::send_dynamics<CartesianMaps, YacCartesianDynamics>(
    const CartesianMaps&, const viewh_constgbx, const YacCartesianDynamics&) const;

template void YacComms::receive_dynamics<CartesianMaps, YacCartesianDynamics>(
    const CartesianMaps&, const YacCartesianDynamics&, const viewh_gbx) const;
