/*
 * ----- CLEO -----
 * File: collisionkinetics.hpp
 * Project: superdrops
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 22nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions for kinetic calculations involved in
 * the collisions of two (real) droplets e.g.
 * used in the probability of coalescence or
 * breakup according to  Low and List 1982(a).
 */

#ifndef COLLISIONKINETCS_HPP
#define COLLISIONKINETCS_HPP

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

KOKKOS_INLINE_FUNCTION
double collision_kinetic_energy(const double r1,
                                const double r2,
                                const double terminalv1,
                                const double terminalv2);
/* returns cke, where cke = collision kinetic energy
as formulated in Low and List 1982(a) eqn 3.1 given
radii r1 and r2 and terminal velocities of droplets */

KOKKOS_INLINE_FUNCTION
double surfenergy(const double radius)
/* returns energy due to surface tension of a single
drop, analogous to equation 4.2 of Low and List 1982 */
{
  /* constant required to calculate surface tension
  energy from dimensionless radius using surface
  tension of water = sigma = 7.28e-2 */
  constexpr double surfconst = 4.0 * 7.28e-2 * Kokkos::numbers::pi *
                               dlc::R0 * dlc::R0; // [J/m^-2]

  return surfconst * radius * radius; // = surfe, droplet surface energy
}

KOKKOS_INLINE_FUNCTION
double total_surfenergy(const double r1, const double r2)
/* returns total energy due to surface tension of
pair of droplets with radii r1 and r2 as in
equation 4.2 of Low and List 1982 */
{
  /* constant required to calculate surface tension
  energy from dimensionless radius using surface
  tension of water = sigma = 7.28e-2 */
  constexpr double surfconst = 4.0 * 7.28e-2 * Kokkos::numbers::pi *
                               dlc::R0 * dlc::R0; // [J/m^-2]
  
  const double rsqrdsum = ((r1 * r1) + (r2 * r2)); // * R0sqrd to convert to [m^2]

  return surfconst * rsqrdsum; // = tot_surfe, total surface energy
}

KOKKOS_INLINE_FUNCTION
double coal_surfenergy(const double r1, const double r2);
/* returns surface energy of single spherical equivalent, ie.
coalesced state of two drops, divided by pi as in
equation 4.3 of Low and List 1982 */


/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_INLINE_FUNCTION
double collision_kinetic_energy(const double r1,
                                const double r2,
                                const double terminalv1,
                                const double terminalv2)
/* returns cke, where cke = collision kinetic energy
as formulated in Low and List 1982(a) eqn 3.1 given
radii r1 and r2 and terminal velocities of droplets */
{
  constexpr double R0cubed = dlc::R0 * dlc::R0 * dlc::R0; // convert r^3 to [m^3]
  constexpr double ckeconst = R0cubed * 2.0 / 3.0 * DC::RHO_L *
                              Kokkos::numbers::pi * dlc::W0;

  const double r1cubed(r1 * r1* r1);
  const double r1_r2cubed((r1 / r2) * (r1 / r2) * (r1 / r2));
  const double rratio(r1cubed / 1 + r1_r2cubed); // * R0cubed to convert to [m^3]

  const double vdiff = terminalv1 - terminalv2; // * dlc::W0 to convert to [m/s]
  const double cke = ckeconst * rratio * vdiff * vdiff;

  return cke;
}


KOKKOS_INLINE_FUNCTION
double coal_surfenergy(const double r1, const double r2)
/* returns surface energy of single spherical equivalent, ie.
coalesced state of two drops, divided by pi as in
equation 4.3 of Low and List 1982 */
{
  /* constant required to calculate surface tension
  energy from dimensionless radius using surface
  tension of water = sigma = 7.28e-2 */
  constexpr double surfconst = 4.0 * 7.28e-2 * Kokkos::numbers::pi *
                               dlc::R0 * dlc::R0; // [J/m^-2]

  const double r1cubed(r1 * r1 * r1);
  const double r2cubed(r2 * r2 * r2);
  const double rcubedsum = r1cubed  + r2cubed;

  const double equiv_surfe = surfconst *
                             Kokkos::pow(rcubedsum, 2.0 / 3.0);

  return equiv_surfe; // coalesced (spherical equivalent) surface energy
}

#endif // COLLISIONKINETCS_HPP
