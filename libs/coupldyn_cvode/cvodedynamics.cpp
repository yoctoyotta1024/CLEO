/*
 * ----- CLEO -----
 * File: cvodedynamics.cpp
 * Project: coupldyn_cvode
 * Created Date: Friday 13th October 2023
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
 * functionality for coupleddyanmics concept for
 * dynamics solver in CLEO where coupling is
 * two-way to cvode adiabatic parcel ODE solver
 */

#include "./cvodedynamics.hpp"

void CvodeDynamics::prepare_to_timestep() const
{
  print_init_ODEdata(step2dimlesstime(couplstep),
                     step2dimlesstime(t_end));
}

void CvodeDynamics::run_dynamics(const unsigned int t_mdl,
                                    const unsigned int t_next) const
{

  previousstates = //TODO
  
  cvode.run_cvodestep(t_mdl, couplstep,
                        step2dimlesstime(t_mdl + couplstep));

  std::cout << "from cvode dyn @ t=" << t_mdl << "\n";
}
