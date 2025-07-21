/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: terminalvelocity.hpp
 * Project: superdrops
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for terminal velocity formulas (for example used by some types of super-droplet
 * motion and collision kernels). Formulas are contained in structures which satisfy the
 * constraints of the VelocityFormula concept
 */

#ifndef LIBS_SUPERDROPS_TERMINALVELOCITY_HPP_
#define LIBS_SUPERDROPS_TERMINALVELOCITY_HPP_

#include <concepts>

#include "../cleoconstants.hpp"
#include "superdrop.hpp"

namespace dlc = dimless_constants;

/**
 * @brief Concept representing a velocity formula for calculating a droplet's terminal velocity.
 *
 * Objects that are of type 'VelocityFormula' take a superdroplet and return something
 * convertible to a double (should be its terminal velocity!).
 */
template <typename V>
concept VelocityFormula = requires(V v, const Superdrop &drop) {
  { v(drop) } -> std::convertible_to<double>;
};

/**
 * @brief Null terminal velocity formula returning zero velocity.
 */
struct NullTerminalVelocity {
  /**
   * @brief Returns zero as the terminal velocity of a droplet.
   *
   * @param drop The super-droplet.
   * @return 0.0
   */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop) const { return 0.0; }
};

/**
 * @brief Terminal velocity formula as in Simmel et al. (2002).
 */
struct SimmelTerminalVelocity {
  /**
   * @brief Returns the mass of a droplet as if it's all water [gramms] to use as 'x' according to
   * Simmel et al. 2002 equation (14).
   *
   * @param radius The radius of the droplet.
   * @return The mass of the droplet as if it's all water [gramms].
   */
  KOKKOS_FUNCTION
  double watermass(const double radius) const;

  /**
   * @brief Returns the (dimensionless) terminal velocity of a droplet according to
   * Simmel et al. (2002).
   *
   * Simmel et al. 2002 formula is a semi-empirical adapted from the work of Gunn and Kinzer (1949)
   * and Beard (1976) and used in Simmel's parmeterisation of Long 1974's hydrodynamic collision
   * kernel. For drops with radius >= 1.74mm, terminal velocity is 9.17m/s.
   *
   * _Note:_ Improvement could be made by following Arabas et al. 2015 and Morrison et al. 2005 in
   * multiplying the terminal velocity by the density ratio, rho_dry0/rho_dry, of dry air under
   * standard conditions (rho_dry0) and in current state (rho_dry).
   *
   * @param drop The superdroplet.
   * @return The (dimensionless) terminal velocity.
   */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop) const;
};

/**
 * @brief Terminal velocity formula based on Rogers and Yau (1989) textbook.
 */
struct RogersYauTerminalVelocity {
  /**
   * @brief Returns the terminal velocity of a superdroplet according to Rogers and Yau (1989).
   *
   * Formula from Rogers and Yau 1989 textbook "a short course in cloud physics" chapter 8. For
   * small droplets formula parameterises Stokes' terminal velocity (valid at low Reynolds numbers
   * for spherical droplets). For drops with radius >= 2mm, terminal velocity is 9m/s.
   *
   * @param drop The superdroplet.
   * @return The (dimensionless) terminal velocity.
   */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop) const;
};

/**
 * @brief Terminal velocity formula based on Rogers et al. (1993).
 */
struct RogersGKTerminalVelocity {
  /**
   * @brief Returns the terminal velocity of a droplet according to Rogers et al. (1993).
   *
   * See "Comparison of Raindrop Size Distributions Measured by Radar Wind Profiler and by Airplane"
   * by  R. R. Rogers, D. Baumgardner, S. A. Ethier, D. A. Carter, and W. L. Ecklund (1993).
   * Formulation is approximation of Gunn and Kinzer (1949) tabulated values.
   *
   * @param drop The superdroplet.
   * @return The (dimensionless) terminal velocity.
   */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop) const;
};

#endif  // LIBS_SUPERDROPS_TERMINALVELOCITY_HPP_
