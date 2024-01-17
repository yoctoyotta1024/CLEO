/*
 * ----- CLEO -----
 * File: breakup_nfrags.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 17th January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * concept and structures for calculating
 * the number of fragments produce.
 * Used e.g. by the DoBreakup struct.
*/

#ifndef BREAKUP_NFRAGS_HPP
#define BREAKUP_NFRAGS_HPP

#include <concepts>

#include <Kokkos_Core.hpp>

#include "./collisionkinetics.hpp"

template <typename F>
concept NFragments = requires(F f,
                              const Superdrop &d1,
                              const Superdrop &d2)
/* Objects that are of type 'NFragments'
take a pair of superdroplets and returns
something convertible to a double (such as
the number of fragments from a breakup event) */
{
  {
    f(d1, d2)
  } -> std::convertible_to<double>;
};

struct ConstNFrags
/* operator always returns constant
number of fragments 'nfrags'. Struct
obeys NFragments concept  */
{
private:
  const double nfrags;
public:
  ConstNFrags(const double nfrags) : nfrags(nfrags) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &d1,
                    const Superdrop &d2) const
  /* always returns constant number of fragments 'nfrags' */
  {
    return nfrags;
  }
};

struct CollisionKineticEnergyNFrags
/* operator returns number of fragments based on collision
kinetic energy. Struct obeys NFragments concept  */
{
public:
  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2) const
  /* returns number of fragments 'nfrags' based on collision
  kinetic energy of droplets according to parameterisation of total
  number of outcomes from Schlottke et al. 2010 (figure 13). 
  Note: nfrags diverges at cke = alpha^(1/beta), so here cke is
  capped at <= ckemax which is value less than alpha^(1/beta) 
  such that nfrags <= 25 */
  {
    constexpr double alpha = 1.5;
    constexpr double beta = 0.135;
    constexpr double ckemax = 16.49789599;

    const auto terminalv = RogersGKTerminalVelocity{};
    const auto cke = collision_kinetic_energy(drop1.get_radius(),
                                              drop2.get_radius(),
                                              terminalv(drop1),
                                              terminalv(drop2));

    const auto cke_cap = double{Kokkos::fmin(ckemax, cke)};                      // limit cke to less than ckemax
    const auto nfrags = double{1.0 / (alpha - Kokkos::pow(cke_cap, beta))};

    return nfrags;
  }
};

#endif // BREAKUP_NFRAGS_HPP
