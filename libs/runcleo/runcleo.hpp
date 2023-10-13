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
#include "./runtimestats.hpp"

unsigned unsigned int next_stepsize(const unsigned int t_mdl,
                                    const CLEOSDM &sdm);

inline unsigned int start_step(const unsigned int t_mdl,
                               const CLEOSDM &sdm,
                               const CoupledDynamics &coupldyn,
                               Gridboxes &gbxs);

inline unsigned int proceed_to_next_step(unsigned int t_mdl,
                                         unsigned int stepsize,
                                         const CLEOSDM &sdm,
                                         const CoupledDynamics &coupldyn,
                                         Gridboxes &gbxs);

inline int timestep_cleo(const unsigned int t_end,
                             const CLEOSDM &sdm,
                             const CoupledDynamics &coupldyn,
                             RunStats &stats,
                             Gridboxes &gbxs,
                             Superdrops &supers);

int run_cleo(const unsigned int t_end,
             const CLEOSDM &sdm,
             const CoupledDynamics &coupldyn)
/* create gridboxes and superdrops, then
timestep CLEO until t_end and with option
to record some runtime statistics */
{
  // generate runtime objects
  RunStats stats;
  Gridboxes gbxs(sdm.generate_gridboxes());
  Superdrops supers(sdm.generate_superdrops());

  // prepare CLEO for timestepping
  coupldyn.prepare_to_timestep();
  sdm.prepare_to_timestep(gbxs, supers);
  stats.pre_timestepping();

  // do timestepping of CLEO from t=0 to t=t_end
  timestep_cleo(t_end, sdm, coupldyn, stats, gbxs, supers);
  stats.post_timestepping();

  // summary of runtime statistics
  stats.summary();

  return 0;
}

inline int timestep_cleo(const unsigned int t_end,
                             const CLEOSDM &sdm,
                             const CoupledDynamics &coupldyn,
                             RunStats &stats,
                             Gridboxes &gbxs,
                             Superdrops &supers)
/* timestep CLEO from t=0 to t=t_end */
{
  unsigned int t_mdl(0);
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    const unsigned int stepsize(start_step(t_mdl, sdm, coupldyn, gbxs));

    /* advance SDM (optionally concurrent to dynamics solver) */
    sdm.run_step(t_mdl, stepsize);

    /* advance dynamics solver (optionally concurrent to SDM) */
    coupldyn.run_step(t_mdl, stepsize);

    /* proceed to next step (in general involves coupling) */
    t_mdl = proceed_to_next_step(t_mdl, stepsize);
  }

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
  if (t_mdl % sdm.couplstep == 0)
  {
    sdm.receive_dynamics(coupldyn, gbxs);
  }

  sdm.obs.observe_startstep(t_mdl, gbxs);

  return next_stepsize(t_mdl, sdm);
}

inline unsigned int proceed_to_next_step(const unsigned int t_mdl,
                                         const unsigned int stepsize,
                                         const CLEOSDM &sdm,
                                         const CoupledDynamics &coupldyn,
                                         Gridboxes &gbxs)
/* returns incremented timestep 't_mdl' of model
by 'stepsize'. Point where communication from
CLEO SDM Gridbox States to coupled dynamics solver
may occur */
{
  if (t_mdl % sdm.couplstep == 0)
  {
    sdm.send_dynamics(coupldyn, gbxs);
  }

  return t_mdl + stepsize;
}

#endif // RUNCLEO_HPP