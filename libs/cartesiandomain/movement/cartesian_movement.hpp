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
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Movement of a superdroplet in a cartesian domain, optionally distributed across
 * more than one MPI process, using certain boundary conditions and a method for superdroplet motion
 */

#ifndef LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOVEMENT_HPP_
#define LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOVEMENT_HPP_

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/movement/cartesian_transport_across_domain.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "gridboxes/transport_across_domain.hpp"
#include "superdrops/motion.hpp"

template <Motion<CartesianMaps> M, BoundaryConditions<CartesianMaps> BCs>
inline auto cartesian_movement(const CartesianMaps &gbxmaps, const M motion,
                               const BCs boundary_conditions) {
  const TransportAcrossDomain<CartesianMaps> auto transport_across_domain =
      CartesianTransportAcrossDomain{};
  return MoveSupersInDomain<CartesianMaps, M, CartesianTransportAcrossDomain, BCs>(
      motion, transport_across_domain, boundary_conditions);
}

#endif  // LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOVEMENT_HPP_
