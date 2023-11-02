/*
 * ----- CLEO -----
 * File: createcartesianmaps.hpp
 * Project: cartesiandomain
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 2nd November 2023
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

#include "../cleoconstants.hpp"
#include "./cartesianmaps.hpp"
#include "initialise/gbxbounds_frombinary.hpp"

CartesianMaps create_cartesian_maps(const unsigned int nspacedims,
                                    std::string_view grid_filename);

#endif // CREATECARTESIANMAPS_HPP