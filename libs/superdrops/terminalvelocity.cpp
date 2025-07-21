/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: terminalvelocity.cpp
 * Project: superdrops
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * some implementation for terminal velocity
 * formulas (used by some types of superdroplet
 * motion and collision kernels). Formulas are
 * contained in structures which satisfy
 * requirements of the VelocityFormula concept
 */

#include "terminalvelocity.hpp"

/**
 * @brief Returns the mass of a droplet as if it's all water [gramms] to use as 'x' according to
 * Simmel et al. 2002 equation (14).
 *
 * @param radius The radius of the droplet.
 * @return The mass of the droplet as if it's all water [gramms].
 */
KOKKOS_FUNCTION
double SimmelTerminalVelocity::watermass(const double radius) const {
  // constant = 4.0/3.0 * pi * density of liquid water
  constexpr double massconst(4.0 / 3.0 * Kokkos::numbers::pi * dlc::Rho_l);

  const auto mass = massconst * radius * radius * radius;

  return mass * dlc::MASS0grams;  // convert dimensionless mass into grams [g]
}

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
double SimmelTerminalVelocity::operator()(const Superdrop &drop) const {
  /* dimensionless values for radii thresholds, For reference, see table 2 of Simmel et al. 2002 */
  constexpr double r1 = 6.7215e-5 / dlc::R0;
  constexpr double r2 = 7.5582e-4 / dlc::R0;
  constexpr double r3 = 1.73892e-3 / dlc::R0;

  /* alpha constants converted from [g^-beta m s^-1] into [g^-beta] units */
  constexpr double VELCONST = (100.0 * dlc::W0);  // convert from [cm/s] into dimensionless velocity
  constexpr double a1 = 457950 / VELCONST;
  constexpr double a2 = 4962 / VELCONST;
  constexpr double a3 = 1732 / VELCONST;
  constexpr double a4 = 917 / VELCONST;

  const auto radius = drop.get_radius();  // dimensionless droplet radius
  if (radius >= r3) {
    return a4;
  } else {
    const auto MASS = watermass(radius);  // droplet mass in grams [g]

    if (radius >= r2) {
      return a3 * Kokkos::pow(MASS, 1.0 / 6.0);  // dimensionless terminal velocity
    } else if (radius >= r1) {
      return a2 * Kokkos::pow(MASS, 1.0 / 3.0);  // dimensionless terminal velocity
    } else {                                     // radius < r1
      return a1 * Kokkos::pow(MASS, 2.0 / 3.0);  // dimensionless terminal velocity
    }
  }
}

/**
 * @brief Returns the terminal velocity of a droplet according to Rogers et al. (1993).
 *
 * See "Comparison of Raindrop Size Distributions Measured by Radar Wind Profiler and by Airplane"
 * by  R. R. Rogers, D. Baumgardner, S. A. Ethier, D. A. Carter, and W. L. Ecklund (1993).
 * Formulation is approximation of Gunn and Kinzer (1949) tabulated values.
 *
 */
KOKKOS_FUNCTION
double RogersYauTerminalVelocity::operator()(const Superdrop &drop) const {
  constexpr double r1 = 3e-5 / dlc::R0;
  constexpr double r2 = 6e-4 / dlc::R0;
  constexpr double r3 = 2e-3 / dlc::R0;

  constexpr double k1 =
      1.19e8 * dlc::R0 * dlc::R0 / dlc::W0;        // k1 in eqn (8.5) converted to [m^-2]
  constexpr double k2 = 8000 * dlc::R0 / dlc::W0;  // k2 in eqn (8.8) converted to [m^-1]
  constexpr double k3 = 201 / dlc::W0;             // k3 in eqn (8.6) in [m^(-1/2)]
  constexpr double k4 = 9 / dlc::W0;               // k4 is max fall speed [dimensionless]

  const auto radius = drop.get_radius();
  if (radius < r1) {
    return k1 * Kokkos::pow(radius, 2.0);  // eqn (8.5)
  } else if (radius < r2) {
    return k2 * radius;  // eqn (8.8)
  } else if (radius < r3) {
    return k3 * Kokkos::pow((radius * dlc::R0), 0.5);  // eqn (8.6)
  } else {                                             // radius >= r3
    return k4;                                         // see text between eqn (8.7) and (8.8)
  }
}

/**
 * @brief Returns the terminal velocity of a droplet according to Rogers et al. (1993).
 *
 * See "Comparison of Raindrop Size Distributions Measured by Radar Wind Profiler and by Airplane"
 * by  R. R. Rogers, D. Baumgardner, S. A. Ethier, D. A. Carter, and W. L. Ecklund (1993).
 * Formulation is approximation of Gunn and Kinzer (1949) tabulated values.
 *
 */
KOKKOS_FUNCTION
double RogersGKTerminalVelocity::operator()(const Superdrop &drop) const {
  constexpr double radius0 =
      3.725 * 1e-4 / dlc::R0;  // dimensionless conversion of D_0 [mm] (to radius)
  constexpr double kcaps =
      2.0 * 4.0 * 1000 * dlc::R0 / dlc::W0;  // dimensionless conversion of K [m/s /mm] (for radius)
  constexpr double smallk =
      -1.0 * 2.0 * 12.0 * 1000 * dlc::R0;    // dimensionless conversion of k [mm^(-1)]
  constexpr double acaps = 9.65 / dlc::W0;   // dimensionless conversion of A [m/s]
  constexpr double bcaps = 10.43 / dlc::W0;  // dimensionless conversion of B [m/s]
  constexpr double ccaps =
      -1.0 * 2.0 * 0.6 * 1000 * dlc::R0;  // dimensionless conversion of C [mm^(-1)]

  const auto radius = drop.get_radius();
  if (radius < radius0) {
    const auto term = double{1.0 - Kokkos::exp(smallk * radius)};
    return term * kcaps * radius;
  } else {
    return acaps - bcaps * Kokkos::exp(ccaps * radius);
  }
}
