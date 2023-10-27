/*
 * ----- CLEO -----
 * File: cvodecomms.hpp
 * Project: coupldyn_cvode
 * Created Date: Friday 27th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 27th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * send and receive dynamics functions for
 * SDM when coupled to the CVODE ODE solver
 */

#ifndef CVODECOMMS_HPP
#define CVODECOMMS_HPP

#include <vector>

#include "../kokkosaliases"
#include "./cvodedynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"

inline void receive_dynamics_from_cvode(const CvodeDynamics &cvode,
                                        const viewh_gbx h_gbxs)
/* update Gridboxes' states using
information received from CVODE dynanmics
solver for  press, temp, qvap and qcond */
{
  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    State &state(h_gbxs(ii).state);
    state.press = cvode.get_pressure(ii);
    state.temp = cvode.get_temperature(ii);
    state.qvap = cvode.get_qvap(ii);
    state.qcond = cvode.get_qcond(ii);
  }
}

void send_dynamics_to_cvode(CvodeDynamics &cvode,
                            const viewh_constgbx h_gbxs);
/* send information from Gridboxes' states
to CVODE dynanmics solver for  temp, qvap
and qcond (excludes press) */

#endif // CVODECOMMS_HPP