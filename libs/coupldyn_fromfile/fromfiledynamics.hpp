/*
 * ----- CLEO -----
 * File: fromfiledynamics.hpp
 * Project: coupldyn_fromfile
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 6th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct obeying coupleddyanmics concept for
 * dynamics solver in CLEO where coupling is
 * one-way and dynamics are read from file
 */

#ifndef FROMFILEDYNAMICS_HPP 
#define FROMFILEDYNAMICS_HPP 

#include <iostream>

#include "initialise/config.hpp"

struct DynamicsVariables
/* contains 1-D vector for each (thermo)dynamic
variable which is ordered by gridbox at every timestep
e.g. press = [p_gbx0(t0), p_gbx1(t0), ,... , p_gbxN(t0), 
p_gbx0(t1), p_gbx1(t1), ..., p_gbxN(t1), ..., p_gbxN(t_end)]
"pos[_X]" gives position of variable in a vector to read
current timestep from for the first gridbox (gbx0)  */
{

}


struct FromFileDynamics
/* type satisfying CoupledDyanmics solver concept
specifically for thermodynamics and wind velocities
that are read from binary files */
{
private:
  const unsigned int interval;

  void run_dynamics(const unsigned int t_mdl) const;
  
public:
  FromFileDynamics(const Config &config,
                   const unsigned int couplstep)
      : interval(couplstep) {}

  auto get_couplstep() const
  {
    return interval;
  }

  void prepare_to_timestep() const;

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void run_step(const unsigned int t_mdl,
                const unsigned int t_next) const
  {
    if (on_step(t_mdl))
    {
      run_dynamics(t_mdl);
    }
  }

};

#endif // FROMFILEDYNAMICS_HPP