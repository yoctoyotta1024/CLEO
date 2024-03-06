/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: collisionkinetics.hpp
 * Project: src
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 6th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions for kinetic calculations involved in
 * the collisions of two (real) droplets e.g.
 * used in the probability of coalescence or
 * breakup according to  Low and List 1982(a).
 */


#ifndef ROUGHPAPER_SRC_COLLISIONKINETICS_HPP_
#define ROUGHPAPER_SRC_COLLISIONKINETICS_HPP_

#include <Kokkos_Core.hpp>

#include "./cleoconstants.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/* returns cke, where cke = collision kinetic energy
as formulated in Low and List 1982(a) eqn 3.1 given
radii r1 and r2 and terminal velocities of droplets */
KOKKOS_FUNCTION double collision_kinetic_energy(const double r1, const double r2,
                                                const double terminalv1, const double terminalv2);

/* returns surface energy of single spherical equivalent, ie.
coalesced state of two drops, divided by pi as in
equation 4.3 of Low and List 1982 */
KOKKOS_FUNCTION double coal_surfenergy(const double r1, const double r2);

/* returns energy due to surface tension of a single
drop, analogous to equation 4.2 of Low and List 1982 */
KOKKOS_INLINE_FUNCTION
double surfenergy(const double radius) {
  const auto rsqrd = double{radius * radius};  // * R0sqrd to convert to [m^2]

  return dlc::surfconst * rsqrd;  // = surfe, droplet surface energy
}

/* returns total energy due to surface tension of
pair of droplets with radii r1 and r2 as in
equation 4.2 of Low and List 1982 */
KOKKOS_INLINE_FUNCTION
double total_surfenergy(const double r1, const double r2) {
  const auto rsqrdsum = double{(r1 * r1) + (r2 * r2)};  // * R0sqrd to convert to [m^2]

  return dlc::surfconst * rsqrdsum;  // = tot_surfe, total surface energy
}

#endif  // ROUGHPAPER_SRC_COLLISIONKINETICS_HPP_
