/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: collisionkinetics.cpp
 * Project: collisions
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * some function implementation for kinetic calculations involved in the collisions of two (real)
 * droplets e.g. used in the probability of coalescence or breakup according to
 * Low and List 1982(a).
 */

#include "./collisionkinetics.hpp"

/**
 * @brief Calculates the collision kinetic energy between two droplets.
 *
 * Returns cke, where cke = collision kinetic energy [Joules] as formulated in Low and
 * List 1982(a) eqn 3.1 given the dimensionless radii, r1 and r2, and the dimensionless
 * terminal velocities of droplets, terminalv1 and terminalv2, respectively.
 *
 * @return The collision kinetic energy [Joules].
 */
KOKKOS_FUNCTION
double collision_kinetic_energy(const double r1, const double r2, const double terminalv1,
                                const double terminalv2) {
  constexpr double R0cubed = dlc::R0 * dlc::R0 * dlc::R0;  // convert r^3 to [m^3]
  constexpr double ckeconst =
      R0cubed * 2.0 / 3.0 * DC::RHO_L * Kokkos::numbers::pi * dlc::W0 * dlc::W0;

  const auto r1cubed = double{r1 * r1 * r1};
  const auto r1_r2cubed = double{(r1 / r2) * (r1 / r2) * (r1 / r2)};
  const auto rratio = double{r1cubed / (1 + r1_r2cubed)};  // * R0cubed to convert to [m^3]

  const auto vdiff = double{terminalv1 - terminalv2};  // * dlc::W0 to convert to [m/s]
  const auto cke = double{ckeconst * rratio * vdiff * vdiff};

  return cke;  // [Joules]
}

/**
 * @brief Calculates the surface tension energy of a coalesced droplet.
 *
 * Returns the surface tension energy of a single spherical droplet, as calculated by equation
 * 4.3 of Low and List 1982, equiavelent to two droplets which coalesce.
 *
 * @param r1 The dimensionless radius of one droplet involed in coalescence.
 * @param r2 The dimensionless radius of the other droplet involed in coalescence.
 * @return The surface tension energy of the superdroplet [Joules].
 */
KOKKOS_FUNCTION
double coal_surfenergy(const double r1, const double r2) {
  const auto r1cubed = double{r1 * r1 * r1};
  const auto r2cubed = double{r2 * r2 * r2};
  const auto rcubedsum = double{r1cubed + r2cubed};

  const auto equiv_surfe = double{dlc::surfconst * Kokkos::pow(rcubedsum, 2.0 / 3.0)};

  return equiv_surfe;  // coalesced surface tension energy [Joules]
}
