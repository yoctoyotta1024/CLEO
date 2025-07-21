/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: lowlistprob.cpp
 * Project: collisions
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for probability of collision event between two (real) droplets according to Low
 * and List . Calculations is contained in structure that satisfies the requirements of the
 * PairProbability concept (see collisions.hpp)
 */

#include "./lowlistprob.hpp"

/* returns coaleff, the coalescence efficiency
of two droplets (given that they have collided)
according to equations (4.5) and (4.6)
Low and List 1982(a) */
KOKKOS_FUNCTION
double LowListCoalProb::coaleff(const Superdrop &drop1, const Superdrop &drop2) const {
  constexpr double aconst = 0.778;
  constexpr double elim = 5e-6;  // total energy limit [J]

  const auto r1 = double{drop1.get_radius()};
  const auto r2 = double{drop2.get_radius()};
  const auto terminalv = RogersGKTerminalVelocity{};

  const auto cke = double{collision_kinetic_energy(r1, r2, terminalv(drop1), terminalv(drop2))};
  const auto surf_t = double{total_surfenergy(r1, r2)};  // [J] S_t
  const auto surf_c = double{coal_surfenergy(r1, r2)};   // [J] S_c
  const auto etot = double{cke + surf_t - surf_c};       // [J] total energy

  if (etot < elim) {
    const auto radiiratio = sizeratio_factor(r1, r2);
    const auto coaleff = double{aconst * radiiratio * expon(etot, surf_c)};

    return coaleff;
  } else {  // coaleff = 0.0
    return 0.0;
  }
}

/* returns the exponential factor in eqn 4.5
Low and List 1982(a) given the total collision energy,
etot [J] and equivalent surface energy, surf_c [J] */
KOKKOS_FUNCTION
double LowListCoalProb::expon(const double etot, const double surf_c) const {
  constexpr double bconst = -2.62e6;  // [J^-2]
  constexpr double sigma = 7.28e-2;   // [J/m^-2]

  const auto exponent = double{bconst * sigma * etot * etot / surf_c};

  return Kokkos::exp(exponent);
}

/* returns factor that takes into account the size
ratio of droplets in eqn 4.5 Low and List 1982(a). */
KOKKOS_FUNCTION
double LowListCoalProb::sizeratio_factor(const double r1, const double r2) const {
  const auto rsmall = double{Kokkos::fmin(r1, r2)};
  const auto rbig = double{Kokkos::fmax(r1, r2)};
  const auto alpha = double{1 + rsmall / rbig};  // alpha = 1 + Ds/Dl

  return 1.0 / (alpha * alpha);  // alpha^(-2)
}
