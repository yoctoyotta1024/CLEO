/*
 * ----- CLEO -----
 * File: cleosdm.cpp
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
 * struct wrapping the core ingredients of the Super-droplet Model
 * (SDM) part of CLEO to enact on super-droplets and gridboxes
 */

#include "./cleosdm.hpp"

void CLEOSDM::prepare_to_timestep(const CoupledDynamics &coupldyn) const
/* prepare CLEO SDM for timestepping */
{
  if (couplstep != coupldyn.get_couplstep())
  {
    const std::string err("coupling timestep of dyanmics "
                          "solver and CLEO SDM are not equal");
    throw std::invalid_argument(err);
  }
}

void CLEOSDM::receive_dynamics(const CoupledDynamics &coupldyn,
                               viewh_gbx h_gbxs) const
/* update Gridboxes' states (on host)
using information received from coupldyn */
{
}

void CLEOSDM::send_dynamics(const CoupledDynamics &coupldyn,
                            viewh_constgbx h_gbxs) const
/* send information from Gridboxes' states
(on host) to coupldyn */
{
}

void CLEOSDM::run_step(const unsigned int t_mdl,
                       const unsigned int stepsize,
                       viewd_gbx d_gbxs,
                       viewd_supers supers) const
/* run CLEO SDM (on device) from time t_mdl to
t_mdl + stepsize with sub-timestepping routine
for super-droplets' movement and microphysics */
{
  unsigned int t_sdm(t_mdl);
  while (t_sdm < t_mdl + stepsize)
  {
    unsigned int t_next(next_sdmstep(t_sdm, stepsize));

    superdrops_movement(t_sdm, d_gbxs, supers);
    sdm_microphysics(t_sdm, t_next, d_gbxs);

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

void CLEOSDM::superdrops_movement(const unsigned int t_sdm,
                                  viewd_gbx d_gbxs,
                                  viewd_supers supers) const
/* move superdroplets (including movement between gridboxes)
according to movesupers struct */
{
  movesupers.run_step(t_sdm, gbxmaps, d_gbxs, supers);
}

void CLEOSDM::sdm_microphysics(const unsigned int t_sdm,
                               const unsigned int t_next,
                               viewd_gbx d_gbxs) const
/* enact SDM microphysics for each gridbox
(using sub-timestepping routine) */
{
  // loop over gbxs for(gbx : gbxs)
  {
    for (unsigned int subt = t_sdm; subt < t_next;
         subt = microphys.next_step(subt))
    {
      microphys.run_step(subt);
    }
  }
}