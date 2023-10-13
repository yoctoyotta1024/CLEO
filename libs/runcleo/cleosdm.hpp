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

struct TimeSteps
{
private:
  unsigned int t_end = 5;
public:
  unsigned int get_t_end() const { return t_end; }
};

struct GridBoxMaps 
{
};

struct Microphys
{
};

struct Motion
{
};

struct Observer
{
};

struct GridBoxes
{
};

struct SuperDrops
{
};

struct CLEOSDM
{
  TimeSteps tsteps; // timesteps for model (e.g. coupling and end time)
  GridBoxMaps gbxmaps; // maps from gridbox indexes to domain coordinates
  Microphys microphys; // microphysical process
  Motion motion; // super-droplets' motion
  Observer obs; // observer

  GridBoxes generate_gridboxes() const;

  SuperDrops generate_superdrops() const;

  int prepare_to_timestep(const GridBoxes &GBxs,
                          const SuperDrops &SDs) const;

  void run_step(const unsigned int t_mdl) const;
};

#endif // CLEOSDM_HPP