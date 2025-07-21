/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesian_motion.hpp
 * Project: movement
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
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

#ifndef LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOTION_HPP_
#define LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOTION_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <cassert>
#include <functional>

#include "../../cleoconstants.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/doubly_periodic_domain.hpp"
#include "gridboxes/predcorrmotion.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

/* wrapper of operator for use of function in PredCorrMotion's CheckBounds type */
struct CartesianCheckBounds {
  /* raise error if superdrop not either out of domain
  or within bounds (ie. lower_bound <= coord < upper_bound) */
  KOKKOS_INLINE_FUNCTION void operator()(const unsigned int idx,
                                         const Kokkos::pair<double, double> bounds,
                                         const double coord) const {
    const bool bad_gbxindex((idx != LIMITVALUES::oob_gbxindex) &&
                            ((coord < bounds.first) | (coord >= bounds.second)));

    assert((!bad_gbxindex) &&
           "SD not in previous gbx nor a neighbour."
           " Try reducing the motion timestep to"
           " satisfy CFL criteria, or use "
           " 'update_ifoutside' to update sd_gbxindex");
  }
};

/* returned type satisfies motion concept for motion of a
superdroplet using a predictor-corrector method to update
a superdroplet's coordinates and then updating it's
sdgbxindex as appropriate for a cartesian domain */
template <VelocityFormula TV>
inline PredCorrMotion<CartesianMaps, TV, CartesianCheckBounds> CartesianMotion(
    const unsigned int motionstep, const std::function<double(unsigned int)> int2time,
    const TV terminalv) {
  return PredCorrMotion<CartesianMaps, TV, CartesianCheckBounds>(motionstep, int2time, terminalv,
                                                                 CartesianCheckBounds{});
}

#endif  // LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_MOTION_HPP_
