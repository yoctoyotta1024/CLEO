/*
 * ----- CLEO -----
 * File: coupleddynamics.cpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for encasing dyanmics solver which coupled to CLEO SDM
 */

#include "./coupleddynamics.hpp"

void CoupledDynamics::prepare_to_timestep() const
{
}

void CoupledDynamics::run_dynamics(const unsigned int t_mdl,
                                   const unsigned int stepsize) const
{
  std::cout << "dyn @ t=" << t_mdl << "\n";
}
