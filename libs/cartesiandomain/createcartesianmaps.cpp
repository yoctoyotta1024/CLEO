/*
 * ----- CLEO -----
 * File: createcartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 1st November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for creating a cartesian maps struct
 * from a GridBoxBoundaries struct containing gridbox's
 * indexes and their coordinate (upper and lower) boundaries
 */

#include "./createmaps_frombinary.hpp"

CartesianMaps create_cartesian_maps(const unsigned int nspacedims,
                                    std::string_view grid_filename)
{
  const GridboxBoundsFromBinary gfb(nspacedims, grid_filename);
  

  return CartesianMaps();
}