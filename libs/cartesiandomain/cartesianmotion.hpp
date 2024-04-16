/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesianmotion.hpp
 * Project: cartesiandomain
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 16th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * the sdgbxindex updated accordingly for a
 * cartesian domain with finite/periodi boundary
 * conditions
 */

#ifndef LIBS_CARTESIANDOMAIN_CARTESIANMOTION_HPP_
#define LIBS_CARTESIANDOMAIN_CARTESIANMOTION_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <cassert>
#include <functional>

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/doubly_periodic_domain.hpp"
#include "gridboxes/predcorrmotion.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

KOKKOS_FUNCTION unsigned int change_if_coord3nghbr(const CartesianMaps &gbxmaps, unsigned int idx,
                                                   Superdrop &drop);

KOKKOS_FUNCTION unsigned int change_if_coord1nghbr(const CartesianMaps &gbxmaps, unsigned int idx,
                                                   Superdrop &drop);

KOKKOS_FUNCTION unsigned int change_if_coord2nghbr(const CartesianMaps &gbxmaps, unsigned int idx,
                                                   Superdrop &drop);

/* wrapper of operator for use of function
in PredCorrMotion's CheckBounds type */
struct CartesianCheckBounds {
  /* raise error if superdrop not either out of domain
  or within bounds (ie. lower_bound <= coord < upper_bound) */
  KOKKOS_INLINE_FUNCTION void operator()(const unsigned int idx,
                                         const Kokkos::pair<double, double> bounds,
                                         const double coord) const {
    const bool bad_gbxindex((idx != outofbounds_gbxindex()) &&
                            ((coord < bounds.first) | (coord >= bounds.second)));

    assert((!bad_gbxindex) &&
           "SD not in previous gbx nor a neighbour."
           " Try reducing the motion timestep to"
           " satisfy CFL criteria, or use "
           " 'update_ifoutside' to update sd_gbxindex");
  }
};

/* wrapper of functions for use in PredCorrMotion's
ChangeToNghbr type for deciding if a superdroplet should move
to a neighbouring gbx in a cartesian domain and then updating the
superdroplet appropriately. Struct has three functions, one
for each direction (coord3 = z, coord1 = x, coord2 = y). For each,
the superdrop's coord is compared to gridbox bounds given by gbxmaps
for the current gbxindex 'idx'. If superdrop coord lies outside
bounds, forward or backward neighbour functions are called to
update sdgbxindex (and possibly other superdrop attributes) */
struct CartesianChangeIfNghbr {
  KOKKOS_INLINE_FUNCTION unsigned int coord3(const CartesianMaps &gbxmaps, unsigned int idx,
                                             Superdrop &drop) const {
    return change_if_coord3nghbr(gbxmaps, idx, drop);
  }

  KOKKOS_INLINE_FUNCTION unsigned int coord1(const CartesianMaps &gbxmaps, unsigned int idx,
                                             Superdrop &drop) const {
    return change_if_coord1nghbr(gbxmaps, idx, drop);
  }

  KOKKOS_INLINE_FUNCTION unsigned int coord2(const CartesianMaps &gbxmaps, unsigned int idx,
                                             Superdrop &drop) const {
    return change_if_coord2nghbr(gbxmaps, idx, drop);
  }
};

/* returned type satisfies motion concept for motion of a
superdroplet using a predictor-corrector method to update
a superdroplet's coordinates and then updating it's
sdgbxindex as appropriate for a cartesian domain */
template <VelocityFormula TV>
inline PredCorrMotion<CartesianMaps, TV, CartesianChangeIfNghbr, CartesianCheckBounds>
CartesianMotion(const unsigned int motionstep, const std::function<double(unsigned int)> int2time,
                const TV terminalv) {
  return PredCorrMotion<CartesianMaps, TV, CartesianChangeIfNghbr, CartesianCheckBounds>(
      motionstep, int2time, terminalv, CartesianChangeIfNghbr{}, CartesianCheckBounds{});
}

#endif  // LIBS_CARTESIANDOMAIN_CARTESIANMOTION_HPP_
