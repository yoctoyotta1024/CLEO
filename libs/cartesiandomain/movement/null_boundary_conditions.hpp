/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: null_boundary_conditions.hpp
 * Project: movement
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 24th Febuary 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * A definition of the Domain Boundary Conditions to use for Cartesian GridBox Maps, Motion of
 * Super-Droplets and MoveSupersInDomain
 */

#ifndef LIBS_CARTESIANDOMAIN_MOVEMENT_NULL_BOUNDARY_CONDITIONS_HPP_
#define LIBS_CARTESIANDOMAIN_MOVEMENT_NULL_BOUNDARY_CONDITIONS_HPP_

#include "../../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "gridboxes/supersindomain.hpp"

// TODO(CB): make boundary conditions a concept and move into gridboxes library
struct NullBoundaryConditions {
  SupersInDomain operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                            const SupersInDomain &allsupers) const {
    return allsupers;
  }
};

#endif  // LIBS_CARTESIANDOMAIN_MOVEMENT_NULL_BOUNDARY_CONDITIONS_HPP_
