/*
 * ----- CLEO -----
 * File: cvodecomms.hpp
 * Project: coupldyn_cvode
 * Created Date: Friday 27th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 29th October 2023
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
#include <array>

#include "../kokkosaliases.hpp"
#include "./cvodedynamics.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/state.hpp"

void receive_dynamics_from_cvode(const CvodeDynamics &cvode,
                                 const viewh_gbx h_gbxs);
/* update Gridboxes' states using
information received from CVODE dynanmics
solver for  press, temp, qvap and qcond */

void send_dynamics_to_cvode(CvodeDynamics &cvode,
                            const viewh_constgbx h_gbxs);
/* send information from Gridboxes' states
to CVODE dynanmics solver for  temp, qvap
and qcond (excludes press) */

#endif // CVODECOMMS_HPP