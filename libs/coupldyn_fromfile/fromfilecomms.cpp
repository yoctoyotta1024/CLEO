/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: fromfilecomms.cpp
 * Project: coupldyn_fromfile
 * Created Date: Tuesday 31st October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * send and receive dynamics functions
 * for SDM when coupled to the fromfile
 * dynamics solver
 */

#include "coupldyn_fromfile/fromfilecomms.hpp"

/* update Gridboxes' states using information received
from FromFileDynamics solver for 1-way coupling to CLEO SDM.
Kokkos::parallel_for([...]) (on host) is equivalent to:
for (size_t ii(0); ii < ngbxs; ++ii){[...]}
when in serial */
template <typename GbxMaps, typename CD>
void FromFileComms::receive_dynamics(const GbxMaps &gbxmaps, const FromFileDynamics &ffdyn,
                                     const viewh_gbx h_gbxs) const {
  const size_t ngbxs(h_gbxs.extent(0));

  Kokkos::parallel_for(
      "receive_dynamics", Kokkos::RangePolicy<HostSpace>(0, ngbxs), [=, *this](const size_t ii) {
        update_gridbox_state(ffdyn, gbxmaps.local_to_global_gridbox_index(ii), h_gbxs(ii));
      });
}

/* updates the state of a gridbox using information
received from FromFileDynamics solver for 1-way
coupling to CLEO SDM */
void FromFileComms::update_gridbox_state(const FromFileDynamics &ffdyn, const size_t ii,
                                         Gridbox &gbx) const {
  State &state(gbx.state);

  state.press = ffdyn.get_press(ii);
  state.temp = ffdyn.get_temp(ii);
  state.qvap = ffdyn.get_qvap(ii);
  state.qcond = ffdyn.get_qcond(ii);

  state.wvel = ffdyn.get_wvel(ii);
  state.uvel = ffdyn.get_uvel(ii);
  state.vvel = ffdyn.get_vvel(ii);
}

template void FromFileComms::send_dynamics<CartesianMaps, FromFileDynamics>(
    const CartesianMaps &, const viewh_constgbx, FromFileDynamics &) const;

template void FromFileComms::receive_dynamics<CartesianMaps, FromFileDynamics>(
    const CartesianMaps &, const FromFileDynamics &, const viewh_gbx) const;
