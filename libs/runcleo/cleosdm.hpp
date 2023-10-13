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

#include <iostream>

#include "./coupleddynamics.hpp"
#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"

struct GridBoxes
{
};

struct SuperDrops
{
};

struct GridBoxMaps 
{
  GridBoxMaps(const Config &config){}
};

struct Microphys
{
  Microphys(const Config &config, const Timesteps &tsteps){}
};

struct Motion
{
  Motion(const Config &config, const Timesteps &tsteps){}
};

struct Observer
{
private:
  unsigned int interval;

public:
  Observer(const Config &config, const Timesteps &tsteps){}
  
  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void observe(const unsigned int t_mdl,
               GridBoxes &gbxs) const
  {
    std::cout << "obs gbxs @ t = " << t_mdl << "\n";
  }

  void observe_startstep(const unsigned int t_mdl,
                         GridBoxes &gbxs) const
  {
    if (on_step(t_mdl))
    {
      observe(t_mdl, gbxs);
    }
  }

  unsigned int get_obsstep() const
  {
    return interval;
  }
};

struct CLEOSDM
{
  GridBoxMaps gbxmaps; // maps from gridbox indexes to domain coordinates
  Microphys microphys; // microphysical process
  Motion motion; // super-droplets' motion
  Observer obs; // observer
  unsigned int couplstep;

  CLEOSDM(const Config &config, const Timesteps &tsteps,
          const unsigned int coupldynstep)
      : gbxmaps(config), microphys(config, tsteps),
        motion(config, tsteps), obs(config, tsteps),
        couplstep(tsteps.get_couplstep())
  {
    if (couplstep != coupldynstep)
    {
      const std::string err("coupling timestep of dyanmics "
                            "solver and CLEO SDM are not equal");
      throw std::invalid_argument(err);
    }
  }

  GridBoxes generate_gridboxes() const;

  SuperDrops generate_superdrops() const;

  int prepare_to_timestep(const GridBoxes &gbxs,
                          const SuperDrops &supers) const;

  void run_step(const unsigned int t_mdl,
                const unsigned int stepsize) const;

  void receive_dynamics(const CoupledDynamics &coupldyn,
                        GridBoxes &gbxs) const;
  
  unsigned int get_couplstep() const { return couplstep; }
};

#endif // CLEOSDM_HPP