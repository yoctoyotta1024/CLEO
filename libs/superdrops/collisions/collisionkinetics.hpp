/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: collisionkinetics.hpp
 * Project: collisions
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions for kinetic calculations involved in the collisions of two (real) droplets e.g.
 * used in the probability of coalescence or breakup according to  Low and List 1982(a).
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_COLLISIONKINETICS_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_COLLISIONKINETICS_HPP_

#include <Kokkos_Core.hpp>

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"
#include "../terminalvelocity.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/**
 * @brief Calculates the collision kinetic energy between two droplets.
 *
 * Returns cke, where cke = collision kinetic energy [Joules] as formulated in Low and
 * List 1982(a) eqn 3.1 given the dimensionless radii, r1 and r2, and the dimensionless
 * terminal velocities of droplets, terminalv1 and terminalv2, respectively.
 *
 * @param r1 The radius of the first superdroplet.
 * @param r2 The radius of the second superdroplet.
 * @param terminalv1 The terminal velocity of the first superdroplet.
 * @param terminalv2 The terminal velocity of the second superdroplet.
 *
 * @return The collision kinetic energy.
 */
KOKKOS_FUNCTION
double collision_kinetic_energy(const double r1, const double r2, const double terminalv1,
                                const double terminalv2);

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
double coal_surfenergy(const double r1, const double r2);

/**
 * @brief Calculates the surface tension energy of a single droplet.
 *
 * Returns the energy due to surface tension of a single droplet, analogous to
 * equation 4.2 of Low and List 1982
 *
 * @param radius The dimensionless radius of the droplet.
 * @return The surface tension energy of the droplet [Joules].
 */
KOKKOS_INLINE_FUNCTION
double surfenergy(const double radius) {
  const auto rsqrd = double{radius * radius};  // * R0sqrd to convert to [m^2]

  return dlc::surfconst * rsqrd;  // = surfe, droplet surface tension energy [Joules].
}

/**
 * @brief Calculates the sum of the surface energy of a pair of droplets.
 *
 * Returns the total energy due to surface tension by summing the individual surface tension
 * energies for a pair of droplets with radii r1 and r2 as in equation 4.2 of Low and List 1982.
 *
 * @param r1 The dimensionless radius of the first superdroplet.
 * @param r2 The dimensionless radius of the second superdroplet.
 * @return The total surface energy of the pair of superdroplets [Joules].
 */
KOKKOS_INLINE_FUNCTION
double total_surfenergy(const double r1, const double r2) {
  const auto rsqrdsum = double{(r1 * r1) + (r2 * r2)};  // * R0sqrd to convert to [m^2]

  return dlc::surfconst * rsqrdsum;  // = tot_surfe, total surface tension energy [Joules].
}

#endif  // LIBS_SUPERDROPS_COLLISIONS_COLLISIONKINETICS_HPP_
