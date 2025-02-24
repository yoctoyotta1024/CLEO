/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: mptrac_movement.hpp
 * Project: mptrac_movement
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
 * Movement of a superdroplet in a cartesian domain using MPTRAC, optionally distributed across
 * more than one MPI process, using certain boundary conditions and a method for superdroplet motion
 */

#ifndef LIBS_CARTESIANDOMAIN_MPTRAC_MOVEMENT_MPTRAC_MOVEMENT_HPP_
#define LIBS_CARTESIANDOMAIN_MPTRAC_MOVEMENT_MPTRAC_MOVEMENT_HPP_

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/mptrac_movement/mptrac_transport_across_domain.hpp"
#include "gridboxes/movesupersindomain.hpp"
#include "superdrops/motion.hpp"

template <Motion<CartesianMaps> M, typename BoundaryConditions>
inline auto mptrac_movement(const CartesianMaps &gbxmaps, const M motion,
                               const BoundaryConditions boundary_conditions) {
  const auto transport_across_domain = MPTRACTransportAcrossDomain{};
  return MoveSupersInDomain<CartesianMaps, M, MPTRACTransportAcrossDomain, BoundaryConditions>(
      motion, transport_across_domain, boundary_conditions);
}

#endif  // LIBS_CARTESIANDOMAIN_MPTRAC_MOVEMENT_MPTRAC_MOVEMENT_HPP_
