/*
 * ----- CLEO -----
 * File: runcleo.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 14th October 2023
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

#include <algorithm>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "../kokkosaliases.hpp"
#include "./cleosdm.hpp"
#include "./coupleddynamics.hpp"
#include "./runtimestats.hpp"

dualview_gbx create_gridboxes();

viewd_supers create_superdrops();

unsigned int next_stepsize(const unsigned int t_mdl,
                           const CLEOSDM &sdm);

inline unsigned int start_step(const unsigned int t_mdl,
                               const CLEOSDM &sdm,
                               const CoupledDynamics &coupldyn,
                               dualview_gbx gbxs);

inline void sdm_step(const unsigned int t_mdl,
                     const unsigned int stepsize,
                     const CLEOSDM &sdm,
                     dualview_gbx gbxs, 
                     viewd_supers supers);

inline void coupldyn_step(const unsigned int t_mdl,
                          const unsigned int stepsize,
                          const CoupledDynamics &coupldyn);

inline unsigned int proceed_to_next_step(unsigned int t_mdl,
                                         unsigned int stepsize,
                                         const CLEOSDM &sdm,
                                         const CoupledDynamics &coupldyn,
                                         dualview_gbx gbxs);

inline int timestep_cleo(const unsigned int t_end,
                             const CLEOSDM &sdm,
                             const CoupledDynamics &coupldyn,
                             RunStats &stats,
                             dualview_gbx gbxs,
                             viewd_supers supers);

int run_cleo(const unsigned int t_end,
             const CLEOSDM &sdm,
             const CoupledDynamics &coupldyn)
/* create gridboxes and superdrops, then
timestep CLEO until t_end and with option
to record some runtime statistics */
{
  // generate runtime objects
  RunStats stats;
  dualview_gbx gbxs(create_gridboxes());
  viewd_supers supers(create_superdrops());

  // prepare CLEO for timestepping
  coupldyn.prepare_to_timestep();
  sdm.prepare_to_timestep(coupldyn);
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
                             dualview_gbx gbxs,
                             viewd_supers supers)
/* timestep CLEO from t=0 to t=t_end */
{
  unsigned int t_mdl(0);
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    const unsigned int stepsize(start_step(t_mdl, sdm, coupldyn, gbxs));

    /* advance SDM (optionally concurrent to dynamics solver) */
    sdm_step(t_mdl, stepsize, sdm, gbxs, supers);
    
    /* advance dynamics solver (optionally concurrent to SDM) */
    coupldyn_step(t_mdl, stepsize, coupldyn);
    
    /* proceed to next step (in general involves coupling) */
    t_mdl = proceed_to_next_step(t_mdl, stepsize, sdm, coupldyn, gbxs);
  }

  return 0;
}

inline unsigned int start_step(const unsigned int t_mdl,
                               const CLEOSDM &sdm,
                               const CoupledDynamics &coupldyn,
                               dualview_gbx gbxs)
/* communication of thermodynamic state from dynamics solver
to CLEO's Gridboxes. Followed by observation. Function then
returns size of step to take given current timestep, t_mdl */
{
  if (t_mdl % sdm.couplstep == 0)
  {
    gbxs.sync_host();
    sdm.receive_dynamics(coupldyn, gbxs.view_host());
    gbxs.modify_host();
  }

  gbxs.sync_host(); 
  sdm.obs.observe_startstep(t_mdl, gbxs.view_host());

  return next_stepsize(t_mdl, sdm);
}

inline void sdm_step(const unsigned int t_mdl,
                     const unsigned int stepsize,
                     const CLEOSDM &sdm,
                     dualview_gbx gbxs, 
                     viewd_supers supers)
/* run CLEO SDM (on device) */
{
  gbxs.sync_device();
  sdm.run_step(t_mdl, stepsize, gbxs.view_device(), supers);
  gbxs.modify_device();
}

inline void coupldyn_step(const unsigned int t_mdl,
                          const unsigned int stepsize,
                          const CoupledDynamics &coupldyn)
/* run coupled dynamics solver (on host) */
{
  coupldyn.run_step(t_mdl, stepsize);
}

inline unsigned int proceed_to_next_step(const unsigned int t_mdl,
                                         const unsigned int stepsize,
                                         const CLEOSDM &sdm,
                                         const CoupledDynamics &coupldyn,
                                         dualview_gbx gbxs)
/* returns incremented timestep 't_mdl' of model
by 'stepsize'. Point where communication from
CLEO SDM Gridbox States to coupled dynamics solver
may occur */
{
  if (t_mdl % sdm.couplstep == 0)
  {
    gbxs.sync_host();
    sdm.send_dynamics(coupldyn, gbxs.view_host());
  }

  return t_mdl + stepsize;
}

#endif // RUNCLEO_HPP