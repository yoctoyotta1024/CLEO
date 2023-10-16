/*
 * ----- CLEO -----
 * File: sdmmethods.hpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 16th October 2023
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

#ifndef SDMMETHODS_HPP
#define SDMMETHODS_HPP

#include <algorithm>
#include <string>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "../kokkosaliases.hpp"
#include "./coupleddynamics.hpp"
#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"
#include "sdmdomain/gridbox.hpp"
#include "sdmdomain/gridboxmaps.hpp"
#include "sdmdomain/movesupersindomain.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/microphysicsprocess.hpp"
#include "observers/observers.hpp"

template <typename CD>
struct SDMMethods
{
private:
  unsigned int couplstep;        // coupled timestep

  unsigned int next_sdmstep(const unsigned int t_sdm,
                            const unsigned int stepsize) const;
  /* given current timestep, t_sdm, work out which event
  (motion or one complete step) is next to occur and return
  the time of the sooner event, (ie. next t_move or t_mdl) */

  void superdrops_movement(const unsigned int t_mdl,
                           viewd_gbx d_gbxs,
                           viewd_supers supers) const;
  /* move superdroplets (including movement between gridboxes)
  according to movesupers struct */

  void sdm_microphysics(const unsigned int t_sdm,
                        const unsigned int t_next,
                        viewd_gbx d_gbxs) const;
  /* enact SDM microphysics for each gridbox
  (using sub-timestepping routine) */

public:
  GridboxMaps gbxmaps;           // maps from gridbox indexes to domain coordinates
  MicrophysicsProcess microphys; // microphysical process
  MoveSupersInDomain movesupers; // super-droplets' motion in domain
  Observer obs;                  // observer

  SDMMethods(const GridboxMaps gbxmaps,
          const MicrophysicsProcess microphys,
          const MoveSupersInDomain movesupers,
          const Observer obs,
          const unsigned int couplstep)
      : couplstep(couplstep),
        gbxmaps(gbxmaps),
        microphys(microphys),
        movesupers(movesupers),
        obs(obs) {}

  unsigned int get_couplstep() const { return couplstep; }
 
  void prepare_to_timestep() const {}
  /* prepare CLEO SDM for timestepping */

  void receive_dynamics(const CD &coupldyn,
                        viewh_gbx h_gbxs) const {}
  /* update Gridboxes' states using information
  received from coupldyn */

  void send_dynamics(const CD &coupldyn,
                     viewh_constgbx h_gbxs) const {}
  /* send information from Gridboxes' states to coupldyn */

  void run_step(const unsigned int t_mdl,
                const unsigned int stepsize,
                viewd_gbx d_gbxs,
                viewd_supers supers) const;
  /* run CLEO SDM (on device) from time t_mdl to
  t_mdl + stepsize with sub-timestepping routine
  for super-droplets' movement and microphysics */
};

template <typename CD>
void SDMMethods<CD>::run_step(const unsigned int t_mdl,
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

template <typename CD>
unsigned int SDMMethods<CD>::next_sdmstep(const unsigned int t_sdm,
                                   const unsigned int stepsize) const
/* given current timestep, t_sdm, work out which event
(motion or one complete step) is next to occur and return
the time of the sooner event, (ie. next t_move or t_mdl) */
{
  const unsigned int next_t_mdl(((t_sdm / stepsize) + 1) * stepsize); // t of next output
  const unsigned int next_t_move(movesupers.next_step(t_sdm));        // t of next sdm movement

  return std::min(next_t_mdl, next_t_move);
}

template <typename CD>
void SDMMethods<CD>::superdrops_movement(const unsigned int t_sdm,
                                  viewd_gbx d_gbxs,
                                  viewd_supers supers) const
/* move superdroplets (including movement between gridboxes)
according to movesupers struct */
{
  movesupers.run_step(t_sdm, gbxmaps, d_gbxs, supers);
}

template <typename CD>
void SDMMethods<CD>::sdm_microphysics(const unsigned int t_sdm,
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

#endif // SDMMETHODS_HPP