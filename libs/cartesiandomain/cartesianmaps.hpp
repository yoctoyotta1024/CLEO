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
#include <Kokkos_UnorderedMap.hpp>

#include "../cleoconstants.hpp"
#include "initialise/config.hpp"
#include "../kokkosaliases.hpp"


// TODO 

namespace dlc = dimless_constants;

struct CartesianMaps 
/* type satisfying GridboxMaps concept specifically
for gridboxes defined on in a cartesian C grid with
equal area and volume for each gridbox */
{
private:
  using kokkosmap = Kokkos::UnorderedMap<unsigned int,
                                         Kokkos::pair<double, double>,
                                         ExecSpace>;
  
  kokkosmap coord3_to_bounds;
  kokkosmap coord1_to_bounds;
  kokkosmap coord2_to_bounds;

public:
  CartesianMaps(const Config &config){}
  /* initilaises coord[X]bounds maps (for X = 1, 2, 3,
  corresponding to x, y, z) to map between gbxindexes and
  gridbox boundaries in a cartiesian domain. The keys of each map
  are the gridbox indexes. The corresponding value of each bounds map
  is that gridbox's {lower boundary, upper boundary}.
  In a non-3D case, boundaries for unused dimensions are the min/max 
  possible (numerical limits), however the area and volume of each
  gridbox remain finite. E.g. In the 0-D case, the maps have 1
  {key, value} for gridbox 0 which are numerical limits, whilst the 
  volume function returns a value determined from the gridfile input */

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  d_coord3bounds(const unsigned int gbxidx) const
  /* returns {lower bound, upper bound}  in coord3
  (z) direction of gridbox with index 'gbxidx'
  on device */
  {
    const auto i(coord3_to_bounds.find(gbxidx)); // index in map of key 'gbxidx'

    return coord3_to_bounds.value_at(i) // value returned by map at index i
  }

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  d_coord1bounds(const unsigned int gbxidx) const
  /* returns {lower bound, upper bound}  in coord1
  (x) direction of gridbox with index 'gbxidx'
  on device */
  {
    const auto i(coord1_to_bounds.find(gbxidx)); // index in map of key 'gbxidx'

    return coord1_to_bounds.value_at(i) // value returned by map at index i
  }


  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<double, double>
  d_coord2bounds(const unsigned int gbxidx) const
  /* returns {lower bound, upper bound}  in coord2
  (y) direction of gridbox with index 'gbxidx'
  on device */
  {
    const auto i(coord2_to_bounds.find(gbxidx)); // index in map of key 'gbxidx'

    return coord2_to_bounds.value_at(i) // value returned by map at index i
  }


  double get_area(const unsigned int gbxidx) const
  {
    return 0.0;
  } // TODO

  double get_volume(const unsigned int gbxidx) const
  {
    return 0.0;
  } // TODO
};

#endif // CARTESIANMAPS_HPP