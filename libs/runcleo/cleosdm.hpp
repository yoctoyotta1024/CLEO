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

struct TimeSteps
{
private:
  unsigned int t_end; // time to end model
  unsigned int couplstep; // stepsize between coupling
public:

  TimeSteps(const Config &config)
      : t_end(5),
        couplstep(2) {}

  unsigned int get_t_end() const { return t_end; }
  unsigned int get_couplstep() const { return couplstep; }
};

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
  Microphys(const Config &config, const TimeSteps &mdlsteps){}
};

struct Motion
{
  Motion(const Config &config, const TimeSteps &mdlsteps){}
};

struct Observer
{
private:
  unsigned int interval;

public:
  Observer(const Config &config, const TimeSteps &mdlsteps){}
  
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

  CLEOSDM(const Config &config, const TimeSteps &mdlsteps)
      : gbxmaps(config), microphys(config, mdlsteps),
        motion(config, mdlsteps), obs(config, mdlsteps) {}

  GridBoxes generate_gridboxes() const;

  SuperDrops generate_superdrops() const;

  int prepare_to_timestep(const GridBoxes &gbxs,
                          const SuperDrops &supers) const;

  void run_step(const unsigned int t_mdl,
                const unsigned int stepsize) const;

  void receive_thermodynamics(const CoupledDynamics &coupldyn,
                              GridBoxes &gbxs) const
  {
  }
};

#endif // CLEOSDM_HPP