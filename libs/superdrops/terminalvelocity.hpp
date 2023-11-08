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
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop) const;
  /* returns (dimensionless) terminal velocity of a superdroplet
  according to Simmel et al. 2002. This is semi-empirical formula
  adapted  from work of Gunn and Kinzer, 1949 and Beard, 1976.
  Used in Simmel's form for Long 1974's hydrodynamic kernel.
  Note: Improvement can be made by following Arabas et al. 2015 and
  Morrison et al. 2005 in multiplying this terminal velocity by the
  density ratio, rho_dry0/rho_dry, of dry air under standard
  conditions (rho_dry0) and in current state (rho_dry). */
};

struct RogersYauTerminalVelocity
{
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop) const;
  /* returns (dimensionless) terminal velocity of a superdroplet
  according to formulas based off Stokes' terminal velocity.
  See Rogers and Yau 1989 textbook "a short course in cloud physics"
  chapter 8. Formula valid at low Reynolds No.s for spherical droplets
  but here formula is used beyond validity. For drops with
  radius >= 2mm, terminal velocity is that of a 2mm
  sized droplet = 9m/s. */
};

#endif // TERMINALVELOCITY_HPP
