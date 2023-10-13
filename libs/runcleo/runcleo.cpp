/*
 * ----- CLEO -----
 * File: runcleo.cpp
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
 * functionality related to timestepping CLEO coupled model
 * (CLEO SDM coupled one-way/two-ways to a Dynamics Solver)
 */


#include "./runcleo.hpp"

unsigned int next_stepsize(const unsigned int t_mdl,
                                    const CLEOSDM &sdm)
/* returns size of next step of model ('onestep')
given current time t_mdl, so that next time
(t_next = t_mdl + onestep) is time of obs or coupl */
{
  const unsigned int couplstep(sdm.get_couplstep());
  const unsigned int obsstep(sdm.obs.get_obsstep());

  const auto next_step = [t_mdl](const unsigned int interval)
  {
    return ((t_mdl / interval) + 1) * interval;
  };

  /* t_next is smaller out of time of next coupl and obs */
  const unsigned int next_coupl(next_step(couplstep));
  const unsigned int next_obs(next_step(obsstep));

  return std::min(next_coupl, next_obs) - t_mdl;
}

