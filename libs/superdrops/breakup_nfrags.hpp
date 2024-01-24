/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: breakup_nfrags.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 24th January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept and structures for calculating
 * the number of fragments produce.
 * Used e.g. by the DoBreakup struct.
 */

#ifndef LIBS_SUPERDROPS_BREAKUP_NFRAGS_HPP_
#define LIBS_SUPERDROPS_BREAKUP_NFRAGS_HPP_

#include <concepts>
#include <Kokkos_Core.hpp>

#include "./collisionkinetics.hpp"

/* Objects that are of type 'NFragments'
take a pair of superdroplets and returns
something convertible to a double (such as
the number of fragments from a breakup event) */
template <typename F>
concept NFragments = requires(F f, const Superdrop &d1, const Superdrop &d2) {
  { f(d1, d2) } -> std::convertible_to<double>;
};

/* operator always returns constant
number of fragments 'nfrags'. Struct
obeys NFragments concept  */
struct ConstNFrags {
 private:
  const double nfrags;

 public:
  explicit ConstNFrags(const double nfrags) : nfrags(nfrags) {}

  /* always returns constant number of fragments 'nfrags' */
  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &d1, const Superdrop &d2) const { return nfrags; }
};

/* operator returns number of fragments based on collision
kinetic energy. Struct obeys NFragments concept  */
struct CollisionKineticEnergyNFrags {
 public:
  /* returns number of fragments 'nfrags' based on collision
  kinetic energy of droplets according to parameterisation of total
  number of outcomes from Schlottke et al. 2010 (figure 13) using
  collision kinetic energy in micro-Joules, with
  modifications:
  1) nfrags diverges at cke = alpha^(1/beta)*1e-6 [Joules], so here
  cke is capped at <= ckemax which is value less than alpha^(1/beta)
  such that nfrags <= 25
  2) as cke -> 0, original formulation tends to
  nfrags -> 2/3 < 2.5 formula multiplied by 3.75
  to make limit nfrags -> 2.5 instead */
  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1, const Superdrop &drop2) const {
    constexpr double alpha = 1.5;
    constexpr double beta = 0.135;
    constexpr double ckemax = 16.49789599/1000000; // [J]

    const auto terminalv = RogersGKTerminalVelocity{};
    const auto cke = collision_kinetic_energy(drop1.get_radius(), drop2.get_radius(),
                                              terminalv(drop1), terminalv(drop2));

    const auto cke_cap = double{Kokkos::fmin(ckemax, cke)};  // limit cke to less than ckemax

    const auto gamma = double{Kokkos::pow(cke_cap*1e6, beta)}; // cke in micro-Joules
    const auto nfrags = double{3.75 / (alpha - gamma)};

    return nfrags;
  }
};

#endif  // LIBS_SUPERDROPS_BREAKUP_NFRAGS_HPP_
