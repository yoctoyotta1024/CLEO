/*
 * ----- CLEO -----
 * File: terminalvelocity.cpp
 * Project: superdrops
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 8th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * implementation of for terminal velocity
 * formulas (used by some types of superdroplet
 * motion and collision kernels). Formulas are
 * contained in structures which satisfy
 * requirements of the VelocityFormula concept
 */

#include "./terminalvelocity.hpp"

KOKKOS_FUNCTION
double SimmelTerminalVelocity::
operator()(const Superdrop &drop) const
/* returns (dimensionless) terminal velocity of a superdroplet
according to Simmel et al. 2002. This is semi-empirical formula
adapted  from work of Gunn and Kinzer, 1949 and Beard, 1976.
Used in Simmel's form for Long 1974's hydrodynamic kernel.
Note: Improvement can be made by following Arabas et al. 2015 and
Morrison et al. 2005 in multiplying this terminal velocity by the
density ratio, rho_dry0/rho_dry, of dry air under standard
conditions (rho_dry0) and in current state (rho_dry). */
{
  /* For reference, see table 2 of Simmel et al. 2002 */
  /* dimensionless values for radii thresholds*/
  constexpr double r1 = 6.7215e-5 / dlc::R0;
  constexpr double r2 = 7.5582e-4 / dlc::R0;
  constexpr double r3 = 1.73892e-3 / dlc::R0;

  /* alpha constants converted from [g^-beta m s^-1] into [g^-beta] units */
  constexpr double MASSCONST = dlc::MASS0grams;  // convert dimensionless mass into grams [g]
  constexpr double VELCONST = (100.0 * dlc::W0); // convert from [cm/s] into dimensionless velocity []
  constexpr double a1 = 457950 / VELCONST;
  constexpr double a2 = 4962 / VELCONST;
  constexpr double a3 = 1732 / VELCONST;
  constexpr double a4 = 917 / VELCONST;

  const double radius(drop.get_radius()); // dimensionless droplet radius []
  const double MASS(drop.mass() * MASSCONST); // droplet mass in grams [g]
  if (radius < r1)
  {
    return a1 * Kokkos::pow(MASS, 2.0 / 3.0); // dimensionless terminal velocity
  }
  else if (radius < r2)
  {
    return a2 * Kokkos::pow(MASS, 1.0 / 3.0); // dimensionless terminal velocity
  }
  else if (radius < r3)
  {
    return a3 * Kokkos::pow(MASS, 1.0 / 6.0); // dimensionless terminal velocity
  }
  else
  {
    return a4;
  }
}

KOKKOS_FUNCTION
double RogersYauTerminalVelocity::
operator()(const Superdrop &drop) const
/* returns (dimensionless) terminal velocity of a superdroplet
according to formulas based off Stokes' terminal velocity.
See Rogers and Yau 1989 textbook "a short course in cloud physics"
chapter 8. Formula valid at low Reynolds No.s for spherical droplets
but here formula is used beyond validity. For drops with
radius >= 2mm, terminal velocity is that of a 2mm
sized droplet = 9m/s. */
{
  constexpr double r1 = 3e-5 / dlc::R0;
  constexpr double r2 = 6e-4 / dlc::R0;
  constexpr double r3 = 2e-3 / dlc::R0;

  constexpr double k1 = 1.19e8 * dlc::R0 * dlc::R0 / dlc::W0; // k1 in eqn (8.5) converted to [m^-2]
  constexpr double k2 = 8000 * dlc::R0 / dlc::W0;             // k2 in eqn (8.8) converted to [m^-1]
  constexpr double k3 = 201 / dlc::W0;                        // k3 in eqn (8.6) in [m^(-1/2)]
  constexpr double k4 = 9 / dlc::W0;                          // k4 is max fall speed [dimensionless]

  const double radius(drop.get_radius());
  if (radius < r1)
  {
    return k1 * Kokkos::pow(radius, 2.0); // eqn (8.5)
  }

  else if (radius < r2)
  {
    return k2 * radius; // eqn (8.8)
  }

  else if (radius < r3)
  {
    return k3 * Kokkos::pow((radius * dlc::R0), 0.5); // eqn (8.6)
  }

  return k4; // dimensionless terminal velocity [see text between eqn (8.7) and (8.8)]
}