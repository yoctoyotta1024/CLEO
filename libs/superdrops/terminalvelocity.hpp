/*
 * ----- CLEO -----
 * File: terminalvelocity.hpp
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
 * Header file for terminal velocity
 * formulas (used by some types of superdroplet
 * motion and collision kernels). Formulas are 
 * contained in structures which satisfy
 * requirements of the VelocityFormula concept
*/

#ifndef TERMINALVELOCITY_HPP
#define TERMINALVELOCITY_HPP

#include <concepts>

#include "../cleoconstants.hpp"
#include "./superdrop.hpp"

namespace dlc = dimless_constants;

template <typename V>
concept VelocityFormula = requires(V v, const Superdrop &drop)
/* Objects that are of type 'VelocityFormula'
take a superdroplet and returns something convertible
to a double (hopefully a velocity!) */
{
  {
    v(drop)
    } -> std::convertible_to<double>;
};

struct NullTerminalVelocity
{
  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop) const
  /* returns 0.0 as terminal velocity of a superdroplet */
  {
    return 0.0;
  }
};

struct SimmelTerminalVelocity
{
  double operator()(const Superdrop &drop) const
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

    double terminal_velocity(0.0);
    if (drop.radius < r1)
    {
      terminal_velocity = a1 * Kokkos::pow((drop.mass()*MASSCONST), 2.0 / 3.0);
    }
    else if (drop.radius < r2)
    {
      terminal_velocity = a2 * Kokkos::pow((drop.mass()*MASSCONST), 1.0 / 3.0);
    }
    else if (drop.radius < r3)
    {
      terminal_velocity = a3 * Kokkos::pow((drop.mass()*MASSCONST), 1.0 / 6.0);
    }
    else
    {
      terminal_velocity = a4;
    }

    return terminal_velocity; // dimensionless velocity
  }
};

#endif // TERMINALVELOCITY_HPP
