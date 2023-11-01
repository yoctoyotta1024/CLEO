/*
 * ----- CLEO -----
 * File: cartesianmaps.hpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
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
 * functions related to creating and using maps to convert
 * between a gridbox indexes and domain coordinates for a 
 * cartesian C grid
 */

#ifndef CARTESIANMAPS_HPP 
#define CARTESIANMAPS_HPP 

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "../cleoconstants.hpp"
#include "initialise/config.hpp"
// TODO 

namespace dlc = dimless_constants;

struct CartesianMaps 
/* type satisfying GridboxMaps concept specifically
for gridboxes defined on in a cartesian C grid with
equal area and volume for each gridbox */
{
private:

public:
  CartesianMaps(const Config &config){}
  /* initilaises coord[X]bounds maps (for X = 1, 2, 3,
  corresponding to x, y, z) to map between gbxindexes and
  gridbox boundaries in a cartiesian domain.
  The keys of each map are the gridbox indexes. The
  corresponding value is that gridbox's {upper boundary, lower boundary}.
  In a non-3D case, coordinates of the gridbox boundaries for unused
  dimensions are the min/max possible doubles of computer (numerical
  limits), however the area and volume remain finite. E.g. In the 0-D
  case, the idx2bounds maps have 1 {key, value} for gridbox 0 which
  are the upper and lower numerical limits, whilst the volume is 
  determined by reading the gridfile */


  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  coord3bounds(const unsigned int gbxidx) const
  {
    return {LIMITVALUES::llim, LIMITVALUES::ulim};
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  coord1bounds(const unsigned int gbxidx) const
  {
    return {LIMITVALUES::llim, LIMITVALUES::ulim}; 
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  coord2bounds(const unsigned int gbxidx) const
  {
    return {LIMITVALUES::llim, LIMITVALUES::ulim};
  }

  double get_area() const { return 0.0; } // TODO

  double get_volume() const { return 0.0; } // TODO
};

#endif // CARTESIANMAPS_HPP