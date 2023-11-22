/*
 * ----- CLEO -----
 * File: collisionkinetics.hpp
 * Project: collisionprobs
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

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"
#include "../terminalvelocity.hpp"

template <VelocityFormula TerminalVelocity>
struct CollisionKinetics
/* calculations involved in the kinetics of
a collision between two superdroplets */
{
private:
  TerminalVelocity terminalv;

  /* constant required to calculate surface tension energy from
  dimensionless radius using surface tension of water = sigma = 7.28e-2 */
  const double surfconst{4.0 * 7.28e-2 * std::numbers::pi * dlc::R0 * dlc::R0}; // [J/m^-2]

public:
  CollisionKinetics(TerminalVelocity tv) : terminalv(tv){};

  double collision_kinetic_energy(const Superdrop &drop1,
                                  const Superdrop &drop2) const
  /* returns cke, where cke = collision kinetic energy
  as formulated in Low and List 1982(a) eqn 3.1 */
  {
    constexpr double R0cubed = dlc::R0 * dlc::R0 * dlc::R0; // convert r^3 to [m^3]
    constexpr double ckeconst = R0cubed * 2.0 / 3.0 * DC::RHO_L *
                                std::numbers::pi * dlc::W0;

    const double r1_r2(drop1.radius / drop2.radius);
    const double rratio = std::pow(drop1.radius, 3.0) /
                          (1 + std::pow(r1_r2, 3.0)); // * R0cubed to convert to [m^3]

    const double vdiff = terminalv(drop1) - terminalv(drop2); // * dlc::W0 to convert to [m/s]
    const double cke = ckeconst * rratio * vdiff * vdiff;

    return cke;
  }

  double surfenergy(const Superdrop &drop) const
  /* returns energy due to surface tension of a single
  drop, analogous to equation 4.2 of Low and List 1982 */
  {
    const double rsqrd = drop.radius * drop.radius; // * R0sqrd to convert to [m^2]
    const double tot_surfe = surfconst * rsqrd;

    return tot_surfe; // total surface energy
  }

  double total_surfenergy(const Superdrop &drop1,
                          const Superdrop &drop2) const
  /* returns total energy due to surface tension of pair
  of drops as in equation 4.2 of Low and List 1982 */
  {
    const double r1(drop1.radius);
    const double r2(drop2.radius);
    const double r2sum = (r1 * r1 + r2 * r2); // * R0sqrd to convert to [m^2]

    const double tot_surfe = surfconst * r2sum;

    return tot_surfe; // total surface energy
  }

  double coal_surfenergy(const Superdrop &drop1,
                         const Superdrop &drop2) const
  /* returns surface energy of single spherical equivalent, ie.
  coalesced state of two drops, divided by pi as in
  equation 4.3 of Low and List 1982 */
  {
    const double r1(drop1.radius);
    const double r2(drop2.radius);
    const double r3sum = std::pow(r1, 3.0) + std::pow(r2, 3.0);

    const double equiv_surfe = surfconst * std::pow(r3sum, 2.0 / 3.0);

    return equiv_surfe; // coalesced (spherical equivalent) surface energy
  }
};

#endif // COLLISIONKINETCS_HPP

