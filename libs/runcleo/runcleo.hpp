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
#include "./sdmmethods.hpp"
#include "./coupleddynamics.hpp"
#include "./runtimestats.hpp"

dualview_gbx create_gridboxes();

viewd_supers create_superdrops();

struct RunCLEO
{
  const SDMMethods &sdm;
  const CoupledDynamics &coupldyn;

  inline int prepare_timestepping() const;

  inline void check_coupling() const;

  inline int timestep_cleo(const unsigned int t_end,
                           RunStats &stats,
                           dualview_gbx gbxs,
                           viewd_supers supers) const;

  inline unsigned int start_step(const unsigned int t_mdl,
                                 dualview_gbx gbxs) const;
  
  unsigned int next_stepsize(const unsigned int t_mdl) const;

  inline void sdm_step(const unsigned int t_mdl,
                       unsigned int stepsize,
                       dualview_gbx gbxs,
                       viewd_supers supers) const;

  inline void coupldyn_step(const unsigned int t_mdl,
                            const unsigned int stepsize) const;

  inline unsigned int proceed_to_next_step(unsigned int t_mdl,
                                           unsigned int stepsize,
                                           dualview_gbx gbxs) const;

  int operator()(const unsigned int t_end) const
  /* create gridboxes and superdrops, then
  timestep CLEO until t_end and with option
  to record some runtime statistics */
  {
    // generate runtime objects
    RunStats stats;
    dualview_gbx gbxs(create_gridboxes());
    viewd_supers supers(create_superdrops());

    // prepare CLEO for timestepping
    prepare_timestepping();
    stats.pre_timestepping();

    // do timestepping of CLEO from t=0 to t=t_end
    timestep_cleo(t_end, stats, gbxs, supers);
    stats.post_timestepping();

    // summary of runtime statistics
    stats.summary();

    return 0;
  }
};

inline int RunCLEO::
    prepare_timestepping() const
{
  coupldyn.prepare_to_timestep();
  sdm.prepare_to_timestep();

  check_coupling();

  return 0;
}

inline void RunCLEO::
    check_coupling() const
{
  if (sdm.get_couplstep() != coupldyn.get_couplstep())
  {
    const std::string err("coupling timestep of dyanmics "
                          "solver and CLEO SDM are not equal");
    throw std::invalid_argument(err);
  }
}

inline int RunCLEO::
    timestep_cleo(const unsigned int t_end,
                  RunStats &stats,
                  dualview_gbx gbxs,
                  viewd_supers supers) const
/* timestep CLEO from t=0 to t=t_end */
{
  unsigned int t_mdl(0);
  while (t_mdl <= t_end)
  {
    /* start step (in general involves coupling) */
    const unsigned int stepsize(start_step(t_mdl, gbxs));

    /* advance SDM (optionally concurrent to dynamics solver) */
    sdm_step(t_mdl, stepsize, gbxs, supers);
    
    /* advance dynamics solver (optionally concurrent to SDM) */
    coupldyn_step(t_mdl, stepsize);
    
    /* proceed to next step (in general involves coupling) */
    t_mdl = proceed_to_next_step(t_mdl, stepsize, gbxs);
  }

  return 0;
}

inline unsigned int RunCLEO::
    start_step(const unsigned int t_mdl,
               dualview_gbx gbxs) const
/* communication of thermodynamic state from dynamics solver
to CLEO's Gridboxes. Followed by observation. Function then
returns size of step to take given current timestep, t_mdl */
{
  if (t_mdl % sdm.get_couplstep() == 0)
  {
    gbxs.sync_host();
    sdm.receive_dynamics(coupldyn, gbxs.view_host());
    gbxs.modify_host();
  }

  gbxs.sync_host(); 
  sdm.obs.observe_startstep(t_mdl, gbxs.view_host());

  return next_stepsize(t_mdl);
}

inline void RunCLEO::
    sdm_step(const unsigned int t_mdl,
             const unsigned int stepsize,
             dualview_gbx gbxs,
             viewd_supers supers) const
/* run CLEO SDM (on device) */
{
  gbxs.sync_device();
  sdm.run_step(t_mdl, stepsize, gbxs.view_device(), supers);
  gbxs.modify_device();
}

inline void RunCLEO::
    coupldyn_step(const unsigned int t_mdl,
                  const unsigned int stepsize) const
/* run coupled dynamics solver (on host) */
{
  coupldyn.run_step(t_mdl, stepsize);
}

inline unsigned int RunCLEO::
    proceed_to_next_step(const unsigned int t_mdl,
                         const unsigned int stepsize,
                         dualview_gbx gbxs) const
/* returns incremented timestep 't_mdl' of model
by 'stepsize'. Point where communication from
CLEO SDM Gridbox States to coupled dynamics solver
may occur */
{
  if (t_mdl % sdm.get_couplstep() == 0)
  {
    gbxs.sync_host();
    sdm.send_dynamics(coupldyn, gbxs.view_host());
  }

  return t_mdl + stepsize;
}

#endif // RUNCLEO_HPP