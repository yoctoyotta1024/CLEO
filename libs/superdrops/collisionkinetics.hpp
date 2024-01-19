/*
 * ----- CLEO -----
 * File: collisionkinetics.hpp
 * Project: superdrops
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 15th January 2024
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
  const auto rsqrd = double{radius * radius}; // * R0sqrd to convert to [m^2]

  return dlc::surfconst * rsqrd; // = surfe, droplet surface energy
}

KOKKOS_INLINE_FUNCTION
double total_surfenergy(const double r1, const double r2)
/* returns total energy due to surface tension of
pair of droplets with radii r1 and r2 as in
equation 4.2 of Low and List 1982 */
{
  const auto rsqrdsum = double{(r1 * r1) + (r2 * r2)}; // * R0sqrd to convert to [m^2]

  return dlc::surfconst * rsqrdsum; // = tot_surfe, total surface energy
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
                              Kokkos::numbers::pi * dlc::W0 * dlc::W0;

  const auto r1cubed = double{r1 * r1 * r1};
  const auto r1_r2cubed = double{(r1 / r2) * (r1 / r2) * (r1 / r2)};
  const auto rratio = double{r1cubed / (1 + r1_r2cubed)}; // * R0cubed to convert to [m^3]

  const auto vdiff = double{terminalv1 - terminalv2}; // * dlc::W0 to convert to [m/s]
  const auto cke = double{ckeconst * rratio * vdiff * vdiff};

  return cke;
}

KOKKOS_INLINE_FUNCTION
double coal_surfenergy(const double r1, const double r2)
/* returns surface energy of single spherical equivalent, ie.
coalesced state of two drops, divided by pi as in
equation 4.3 of Low and List 1982 */
{
  const auto r1cubed = double{r1 * r1 * r1};
  const auto r2cubed = double{r2 * r2 * r2};
  const auto rcubedsum = double{r1cubed + r2cubed};

  const auto equiv_surfe = double{dlc::surfconst *
                           Kokkos::pow(rcubedsum, 2.0 / 3.0)};

  return equiv_surfe; // coalesced (spherical equivalent) surface energy
}

#endif // COLLISIONKINETCS_HPP
