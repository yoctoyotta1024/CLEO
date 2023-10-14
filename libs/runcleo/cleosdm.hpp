/*
 * ----- CLEO -----
 * File: cleosdm.hpp
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

#ifndef CLEOSDM_HPP
#define CLEOSDM_HPP

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

struct CLEOSDM
{
private:
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
  unsigned int couplstep;        // coupled timestep
  
  CLEOSDM(const GridboxMaps gbxmaps,
          const MicrophysicsProcess microphys,
          const MoveSupersInDomain movesupers,
          const Observer obs,
          const unsigned int couplstep)
      : gbxmaps(gbxmaps), microphys(microphys),
        movesupers(movesupers), obs(obs),
        couplstep(couplstep) {}

  unsigned int get_couplstep() const { return couplstep; }
 
  void prepare_to_timestep(const CoupledDynamics &coupldyn) const;
  /* prepare CLEO SDM for timestepping */

  void receive_dynamics(const CoupledDynamics &coupldyn,
                        viewh_gbx h_gbxs) const;
  /* update Gridboxes' states using information
  received from coupldyn */

  void send_dynamics(const CoupledDynamics &coupldyn,
                     viewh_constgbx h_gbxs) const;
  /* send information from Gridboxes' states to coupldyn */

  void run_step(const unsigned int t_mdl,
                const unsigned int stepsize,
                viewd_gbx d_gbxs,
                viewd_supers supers) const;
  /* run CLEO SDM (on device) from time t_mdl to
  t_mdl + stepsize with sub-timestepping routine
  for super-droplets' movement and microphysics */
};

#endif // CLEOSDM_HPP