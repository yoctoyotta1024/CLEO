/*
 * ----- CLEO -----
 * File: cartesianmaps.hpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions related to creating and using maps to convert
 * between a gridbox indexes and domain coordinates for a 
 * cartesian C grid
 */

#ifndef CARTESIANMAPS_HPP 
#define CARTESIANMAPS_HPP 

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "initialise/config.hpp"

struct CartesianMaps 
/* type satisfying GridboxMaps concept specifically
for gridboxes defined on in a cartesian C grid */
{
private:

public:
  CartesianMaps(const Config &config){}

  double volume(const unsigned int gbxidx) const
  {
    return 0.0;
  }

  Kokkos::pair<double, double>
  coord3bounds(const unsigned int gbxidx) const
  {
    return {0.0, 1.0};
  }

  Kokkos::pair<double, double>
  coord1bounds(const unsigned int gbxidx) const
  {
    return {0.0, 1.0}; 
  }

  Kokkos::pair<double, double>
  coord2bounds(const unsigned int gbxidx) const
  {
    return {0.0, 1.0};
  }

};

#endif // CARTESIANMAPS_HPP