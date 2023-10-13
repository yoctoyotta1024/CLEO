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

void start_step(const unsigned int t_mdl)
{}

unsigned int next_step(unsigned int t_mdl)
{
  return ++t_mdl;
}

int run_cleo(const CLEOSDM &sdm, const CoupledDynamics &coupldyn)
{
  // generate runtime objects
  RunStats stats;
  GridBoxes GBxs(sdm.generate_gridboxes());
  SuperDrops SDs(sdm.generate_superdrops());

  // prepare CLEO for timestepping
  coupldyn.prepare_to_timestep();
  sdm.prepare_to_timestep(GBxs, SDs);
  stats.pre_timestepping();
  
  // timestep CLEO from t=0 to t=t_end
  unsigned int t_mdl = 0;
  while (t_mdl <= sdm.tsteps.get_t_end())
  {
    /* start step (in general involves coupling) */
    start_step(t_mdl);

    /* advance SDM (optionally concurrent to dynamics solver) */
    sdm.run_step(t_mdl);

    /* advance dynamics solver (optionally concurrent to SDM) */
    coupldyn.run_step(t_mdl);

    /* proceed to next step (in general involves coupling) */
    t_mdl = next_step(t_mdl);
  }
  stats.post_timestepping();
  
  // summary of runtime statistics
  stats.summary();

  return 0;
}

#endif // RUNCLEO_HPP