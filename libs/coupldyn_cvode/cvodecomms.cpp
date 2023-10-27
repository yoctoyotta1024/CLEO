/*
 * ----- CLEO -----
 * File: cvodecomms.cpp
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

#include "./cvodecomms.hpp"

void send_dynamics_to_cvode(CvodeDynamics &cvode,
                            const viewh_constgbx h_gbxs)
/* send information from Gridboxes' states
to CVODE dynanmics solver for  temp, qvap
and qcond (excludes press) */
{
  constexpr size_t NVARS = 4;

  std::vector<double> delta_y(NVARS * ngbxs, 0.0);
  for (size_t ii(0); ii<ngbxs; ++ii)
  {
    const size_t jj(NVARS * ii)
    const State state(h_gbxs(ii).state);

    delta_y.at(jj + 1) = state.temp - previousstates.at(jj + 1);
    delta_y.at(jj + 2) = state.qvap - previousstates.at(jj + 2);
    delta_y.at(jj + 3) = state.qcond - previousstates.at(jj + 3);
  }

  std::vector<double> nodelta(NVARS * ngbxs, 0.0);
  if (delta_y != nodelta)
  {
    cvode.reinitialise(cvode.get_time(), delta_y);
  }
}