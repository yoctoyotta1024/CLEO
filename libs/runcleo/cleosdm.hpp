/*
 * ----- CLEO -----
 * File: cleosdm.hpp
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

#ifndef CLEOSDM_HPP
#define CLEOSDM_HPP

#include <string>
#include <stdexcept>

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
  GridboxMaps gbxmaps; // maps from gridbox indexes to domain coordinates
  MicrophysicsProcess microphys; // microphysical process
  MoveSupersInDomain movesupers; // super-droplets' motion in domain
  Observer obs; // observer
  unsigned int couplstep;

  CLEOSDM(const Config &config, const Timesteps &tsteps,
          const unsigned int coupldynstep)
      : gbxmaps(config), microphys(),
        movesupers(config, tsteps), obs(),
        couplstep(tsteps.get_couplstep())
  {
    if (couplstep != coupldynstep)
    {
      const std::string err("coupling timestep of dyanmics "
                            "solver and CLEO SDM are not equal");
      throw std::invalid_argument(err);
    }
  }

  Gridboxes generate_gridboxes() const;

  Superdrops generate_superdrops() const;

  void prepare_to_timestep(const Gridboxes &gbxs,
                          const Superdrops &supers) const;

  void run_step(const unsigned int t_mdl,
                const unsigned int stepsize) const;

  void receive_dynamics(const CoupledDynamics &coupldyn,
                        Gridboxes &gbxs) const;

  void send_dynamics(const CoupledDynamics &coupldyn,
                     Gridboxes &gbxs) const;

  unsigned int get_couplstep() const { return couplstep; }
};

#endif // CLEOSDM_HPP