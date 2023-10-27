/*
 * ----- CLEO -----
 * File: cvodedynamics.hpp
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
 * struct obeying coupleddyanmics concept for
 * dynamics solver in CLEO where coupling is
 * two-way to cvode adiabatic parcel ODE solver
 */

#ifndef CVODEDYNAMICS_HPP 
#define CVODEDYNAMICS_HPP 

#include <iostream>

#include "initialise/config.hpp"

struct CvodeDynamics
/* type satisfying CoupledDyanmics solver
concept specifically for thermodynamics
of adiabatically expanding parcel (0-D) */
{
private:
  const unsigned int interval;
  std::vector<double> previousstates; // holds states press, temp, qvap and qcond before timestep iterated

  void run_dynamics(const unsigned int t_mdl,
                    const unsigned int t_next) const;
  
public:
  CvodeDynamics(const Config &config,
                const unsigned int couplstep)
      : interval(couplstep)
  {
    /* CVODE thermodynamics solver */
    CvodeThermoSolver cvode(config,
                            initcvodethermo(sdm.gbxmaps.ngridboxes, config));
  }

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
      run_dynamics(t_mdl, t_next);
    }
  }

};

#endif // CVODEDYNAMICS_HPP