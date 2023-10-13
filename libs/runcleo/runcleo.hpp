/*
 * ----- CLEO -----
 * File: runcleo.hpp
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
 * Generic functions for timestepping CLEO coupled model
 * (CLEO SDM coupled one-way/two-ways to a Dynamics Solver)
 */

#ifndef RUNCLEO_HPP 
#define RUNCLEO_HPP 

#include "./cleosdm.hpp"
#include "./coupleddynamics.hpp"
#include "./runstats.hpp"

unsigned int next_stepsize(const unsigned int t_mdl,
                                  const CLEOSDM &sdm);

inline unsigned int start_step(const unsigned int t_mdl,
                               const CLEOSDM &sdm,
                               const CoupledDynamics &coupldyn,
                               Gridboxes &gbxs);

inline unsigned int proceed_to_next_step(unsigned int t_mdl);

int run_cleo(const unsigned int t_end,
             const CLEOSDM &sdm,
             const CoupledDynamics &coupldyn)
{
  // generate runtime objects
  RunStats stats;
  Gridboxes gbxs(sdm.generate_gridboxes());
  Superdrops supers(sdm.generate_superdrops());

  // prepare CLEO for timestepping
  coupldyn.prepare_to_timestep();
  sdm.prepare_to_timestep(gbxs, supers);
  stats.pre_timestepping();
  
  // timestep CLEO from t=0 to t=t_end
  unsigned int t_mdl(0);
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    const unsigned int stepsize = start_step(t_mdl, sdm, coupldyn, gbxs);

    /* advance SDM (optionally concurrent to dynamics solver) */
    sdm.run_step(t_mdl, stepsize);

    /* advance dynamics solver (optionally concurrent to SDM) */
    coupldyn.run_step(t_mdl, stepsize);

    /* proceed to next step (in general involves coupling) */
    t_mdl = proceed_to_next_step(t_mdl);
  }
  stats.post_timestepping();
  
  // summary of runtime statistics
  stats.summary();

  return 0;
}

inline unsigned int start_step(const unsigned int t_mdl,
                               const CLEOSDM &sdm,
                               const CoupledDynamics &coupldyn,
                               Gridboxes &gbxs)
/* communication of thermodynamic state from dynamics solver
to CLEO's Gridboxes. Followed by observation. Function then 
returns size of step to take given current timestep, t_mdl */
{
  sdm.receive_dynamics(coupldyn, gbxs);

  sdm.obs.observe_startstep(t_mdl, gbxs);

  return next_stepsize(t_mdl, sdm);
}

inline unsigned int proceed_to_next_step(unsigned int t_mdl)
{
  return ++t_mdl;
}

#endif // RUNCLEO_HPP