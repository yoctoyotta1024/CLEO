/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: null_boundary_conditions.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 18th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * A definition of the Domain Boundary Conditions to use for Cartesian GridBox Maps, Motion of
 * Super-Droplets and MoveSupersInDomain
 */

#ifndef LIBS_CARTESIANDOMAIN_NULL_BOUNDARY_CONDITIONS_HPP_
#define LIBS_CARTESIANDOMAIN_NULL_BOUNDARY_CONDITIONS_HPP_

#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"

struct NullBoundaryConditions {
  void operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                  const viewd_supers totsupers) const {}
};

#endif  // LIBS_CARTESIANDOMAIN_NULL_BOUNDARY_CONDITIONS_HPP_
