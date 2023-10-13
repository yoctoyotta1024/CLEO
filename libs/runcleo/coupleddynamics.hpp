/*
 * ----- CLEO -----
 * File: coupleddynamics.hpp
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
 * struct for encasing dyanmics solver which coupled to CLEO SDM
 */

#ifndef COUPLEDDYNAMICS_HPP 
#define COUPLEDDYNAMICS_HPP 

#include <iostream>

struct CoupledDynamics
{
  int prepare_to_timestep() const;

  void run_step(const unsigned int t_mdl) const;
};

#endif // COUPLEDDYNAMICS_HPP