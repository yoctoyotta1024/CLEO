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
private:
  const unsigned int interval;

public:
  CoupledDynamics(const unsigned int couplstep)
      : interval(couplstep) {}

  void prepare_to_timestep() const;

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void run_dynamics(const unsigned int t_mdl,
                    const unsigned int stepsize) const;

  void run_step(const unsigned int t_mdl,
                const unsigned int stepsize) const
  {
    if (on_step(t_mdl))
    {
      run_dynamics(t_mdl, stepsize);
    }
  }

  unsigned int get_couplstep() const
  {
    return interval;
  }
};

#endif // COUPLEDDYNAMICS_HPP