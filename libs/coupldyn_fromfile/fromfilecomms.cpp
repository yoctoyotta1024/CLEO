/*
 * ----- CLEO -----
 * File: fromfilecomms.cpp
 * Project: coupldyn_fromfile
 * Created Date: Tuesday 31st October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * send and receive dynamics functions
 * for SDM when coupled to the fromfile
 * dynamics solver
 */

#include "./fromfilecomms.hpp"

template <typename CD>
void FromFileComms::receive_dynamics(const FromFileDynamics &ffdyn,
                                     const viewh_gbx h_gbxs) const
/* update Gridboxes' states using information
received from FromFileDynamics solver for
1-way coupling to CLEO SDM */
{
  const size_t ngbxs(h_gbxs.extent(0));
  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    State &state(h_gbxs(ii).state);

    state.press = ffdyn.get_press(ii);
    state.temp = ffdyn.get_temp(ii);
    state.qvap = ffdyn.get_qvap(ii);
    state.qcond = ffdyn.get_qcond(ii);

    state.wvel = ffdyn.get_wvel(ii);
    state.uvel = ffdyn.get_uvel(ii);
    state.vvel = ffdyn.get_vvel(ii);
  }
}

template void FromFileComms::
    receive_dynamics<FromFileDynamics>(const FromFileDynamics &,
                                       const viewh_gbx) const;