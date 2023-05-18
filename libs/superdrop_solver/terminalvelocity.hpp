// Author: Clara Bayley
// File: terminalvelocity.hpp
/* Header file for terminal velocity
formulas (used in sedimentation method
and some coalescence kernels). Formulas
contained in structures to satisfy
requirements for using VelocityFormula
concept */

#ifndef TERMINALVELOCITY_HPP
#define TERMINALVELOCITY_HPP

#include <concepts>
#include <cmath>

#include "../claras_SDconstants.hpp"
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
  double operator()(const Superdrop &drop) const
  /* returns 0.0 as terminal velocity of a superdroplet */
  {
    return 0.0;
  }
};

struct RogersYauTerminalVelocity
{
  double operator()(const Superdrop &drop) const
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
    constexpr double k4 = 9 / dlc::W0;      // k4 is max fall speed [dimensionless]

    if (drop.radius < r1)
    {
      return k1 * pow(drop.radius, 2.0); // eqn (8.5)
    }

    else if (drop.radius < r2)
    {
      return k2 * drop.radius; // eqn (8.8)
    }

    else if (drop.radius < r3)
    {
      return k3 * pow((drop.radius*dlc::R0), 0.5); // eqn (8.6)
    }

    return k4; // see text between eqn (8.7) and (8.8)
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
    constexpr double MASSCONST = (dlc::R0 * dlc::R0 * dlc::R0) * dlc::RHO0 * 1000; // convert dimensionless mass into grams [g]
    constexpr double VELCONST = (100.0 * dlc::W0);                                 // convert from [cm/s] into dimensionless velocity []
    constexpr double a1 = 457950 / VELCONST;
    constexpr double a2 = 4962 / VELCONST;
    constexpr double a3 = 1732 / VELCONST;
    constexpr double a4 = 917 / VELCONST;

    double terminal_velocity(0.0);
    if (drop.radius < r1)
    {
      terminal_velocity = a1 * pow((drop.mass()*MASSCONST), 2.0 / 3.0);
    }
    else if (drop.radius < r2)
    {
      terminal_velocity = a2 * pow((drop.mass()*MASSCONST), 1.0 / 3.0);
    }
    else if (drop.radius < r3)
    {
      terminal_velocity = a3 * pow((drop.mass()*MASSCONST), 1.0 / 6.0);
    }
    else
    {
      terminal_velocity = a4;
    }

    return terminal_velocity; // from cm/s into dimensionless velocity
  }
};

#endif // TERMINALVELOCITY_HPP