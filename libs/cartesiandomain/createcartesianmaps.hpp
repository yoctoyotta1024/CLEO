/*
 * ----- CLEO -----
 * File: createcartesianmaps.hpp
 * Project: cartesiandomain
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions for creating a cartesian maps struct
 * from a GbxBoundsFromBinary struct containing
 * vectors of gridbox's indexes and their
 * coordinate (upper and lower) boundaries
 */


#ifndef CREATECARTESIANMAPS_HPP
#define CREATECARTESIANMAPS_HPP

#include <string_view>
#include <vector>
#include <stdexcept>

#include "Kokkos_Core.hpp"
#include "Kokkos_Pair.hpp"

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./cartesianmaps.hpp"
#include "initialise/gbxbounds_frombinary.hpp"

/* NOTE: boundary conditions of domain are defined as:
  z: FINITE   (see cartesian_znghbrs, coord3_beyondzdown and coord3_beyondzup) // TODO beyond funcs
  x: PERIODIC (see cartesian_xnghbrs, coord1_beyondxbehind and coord1_beyondxinfront)
  y: PERIODIC (see cartesian_ynghbrs, coord2_beyondyleft and coord2_beyondyright)
*/

CartesianMaps create_cartesian_maps(const unsigned int ngbxs,
                                    const unsigned int nspacedims,
                                    std::string_view grid_filename);
/* creates cartesian maps instance using gridbox bounds read from
gridfile for a 0-D, 1-D, 2-D or 3-D model with periodic or finite 
boundary conditions (see note above). In a non-3D case, boundaries
and neighbours maps for unused dimensions are 'null'
(ie. return numerical limits), however the area and volume of each
gridbox remains finite. E.g. In the 0-D case, the bounds maps all
have 1 {key, value} where key=gbxidx=0 and value = {max, min}
numerical limits, meanwhile volume function returns a value determined
from the gridfile 'grid_filename' */

#endif // CREATECARTESIANMAPS_HPP