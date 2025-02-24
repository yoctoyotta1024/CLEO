/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesian_movement.hpp
 * Project: movement
 * Created Date: Monday 24th Febuary 2025
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
 * Movement of a superdroplet in a cartesian domain, optionally distributed across
 * more than one MPI process, using certain boundary onditions and the predictor corretor
 * method for superdroplet motion
 */

#ifndef LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOVEMENT_HPP_
#define LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOVEMENT_HPP_

#include "cartesiandomain/movement/cartesian_transport_across_domain.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "superdrops/motion.hpp"

template <Motion<CartesianMaps> M, typename BoundaryConditions>
inline auto cartesian_movement(const CartesianMaps &gbxmaps,
  const M motion, const BoundaryConditions boundary_conditions) {
    const auto transport_across_domain = CartesianTransportAcrossDomain{};
  return MoveSupersInDomain<CartesianMaps, M, CartesianTransportAcrossDomain, BoundaryConditions>(
    motion, transport_across_domain, boundary_conditions);
}

#endif // LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOVEMENT_HPP_ 