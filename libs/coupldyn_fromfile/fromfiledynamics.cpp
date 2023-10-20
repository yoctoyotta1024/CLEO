/*
 * ----- CLEO -----
 * File: fromfiledynamics.cpp
 * Project: coupldyn_fromfile
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for dynamics solver in CLEO
 * where coupling is one-way and dynamics
 * are read from file
 */

#include "./fromfiledynamics.hpp"

void FromFileDynamics::prepare_to_timestep() const
{
}

void FromFileDynamics::run_dynamics(const unsigned int t_mdl,
                                    const unsigned int t_next) const
{
  std::cout << "from file dyn @ t=" << t_mdl << "\n";
}
