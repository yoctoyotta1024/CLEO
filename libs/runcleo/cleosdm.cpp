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

GridBoxes CLEOSDM::generate_gridboxes() const
{
  return GridBoxes{};
}

SuperDrops CLEOSDM::generate_superdrops() const
{
  return SuperDrops{};
}

int CLEOSDM::prepare_to_timestep(const GridBoxes &gbxs,
                                 const SuperDrops &supers) const
{
  return 0;
}

void CLEOSDM::receive_dynamics(const CoupledDynamics &coupldyn,
                                     GridBoxes &gbxs) const
{
}

void CLEOSDM::run_step(const unsigned int t_mdl,
                       const unsigned int stepsize) const
{
  std::cout << "SDM Call @ t=" << t_mdl << " (SDM substepping) \n";
}