/*
 * ----- CLEO -----
 * File: cleosdm.cpp
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
 * struct wrapping the core ingredients of the Super-droplet Model
 * (SDM) part of CLEO to enact on super-droplets and gridboxes
 */

#include "./cleosdm.hpp"

void CLEOSDM::prepare_to_timestep(const Gridboxes &gbxs,
                                  const Superdrops &supers) const
/* prepare CLEO SDM for timestepping */
{
}

void CLEOSDM::receive_dynamics(const CoupledDynamics &coupldyn,
                               Gridboxes &gbxs) const
/* update Gridboxes' states using information
received from coupldyn */
{
}

void CLEOSDM::send_dynamics(const CoupledDynamics &coupldyn,
                            Gridboxes &gbxs) const
/* send information from Gridboxes' states to coupldyn */
{
}

void CLEOSDM::run_step(const unsigned int t_mdl,
                       const unsigned int stepsize,
                       Gridboxes &gbxs,
                       Superdrops &supers) const
/* run CLEO SDM from time t_mdl to t_mdl + stepsize with
sub-timestepping routine for super-droplets' movement
and microphysics */
{
  unsigned int t_sdm(t_mdl);
  while (t_sdm < t_mdl + stepsize)
  {
    unsigned int t_next(next_sdmstep(t_sdm, stepsize));

    superdrops_movement(t_sdm, gbxs, supers);
    sdm_microphysics(t_sdm, t_next, gbxs);

    t_sdm = t_next;
  }
}

unsigned int CLEOSDM::next_sdmstep(const unsigned int t_sdm,
                                   const unsigned int stepsize) const
/* given current timestep, t_sdm, work out which event
(motion or one complete step) is next to occur and return
the time of the sooner event, (ie. next t_move or t_mdl) */
{
  const unsigned int next_t_mdl(((t_sdm / stepsize) + 1) * stepsize); // t of next output
  const unsigned int next_t_move(movesupers.next_step(t_sdm));        // t of next sdm movement

  return std::min(next_t_mdl, next_t_move);
}

void CLEOSDM::superdrops_movement(const unsigned int t_mdl,
                                  Gridboxes &gbxs,
                                  Superdrops &supers) const
/* move superdroplets (including movement between gridboxes)
according to movesupers struct */
{
  movesupers.run_step(t_sdm, gbxmaps, gbxs, supers);
}

void CLEOSDM::sdm_microphysics(const unsigned int t_sdm,
                               const unsigned int t_next,
                               Gridboxes &gbxs) const
/* enact SDM microphysics for each gridbox
(using sub-timestepping routine) */
{
  // loop over gbxs for(gbx : gbxs)
  {
    for (unsigned int subt = t_sdm; subt < t_next;
         subt = microphys.next_step(subt))
    {
      microphys.run_step();
    }
  }
}