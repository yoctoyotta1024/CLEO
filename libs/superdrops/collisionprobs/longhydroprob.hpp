/*
 * ----- CLEO -----
 * File: longhydroprob.hpp
 * Project: collisionprobs
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 14th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Probability of collision-coalescence event
 * between two (real) droplets using the
 * hydrodynamic (i.e. gravitational) kernel
 * according to Simmel et al. 2002's formulation
 * of Long's Hydrodynamic Kernel. Probability
 * calculations are contained in structures
 * that satisfy the requirements of the
 * PairProbability concept (see collisions.hpp)
 */

#ifndef LONGHYDROPROB_HPP
#define LONGHYDROPROB_HPP

#include <Kokkos_Core.hpp>

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"
#include "../terminalvelocity.hpp"
#include "./hydrodynamicprob.hpp"

namespace dlc = dimless_constants;

struct LongHydroProb
/* Opertator returns the collision-coalescence probability
  given the efficiency factor, eff = eff(drop1, drop2),
  from Long's hydrodynamic kernel according to Simmel et al. 2002.
*/
{
private:
  const HydrodynamicProb<SimmelTerminalVelocity> hydroprob;
  const double coaleff;

  KOKKOS_INLINE_FUNCTION
  double kerneleff(const Superdrop &drop1,
                   const Superdrop &drop2) const;
  /* returns the efficiency of collision-coalescence, eff, according to
  equations 12 and 13 of Simmel et al. 2002). eff = eff(R,r) where R>r.
  eff = colleff(R,r) * coaleff(R,r). Usually it's assumed that
  coaleff(R,r) = 1, ie. eff = colleff, which also means that for
  collisions where R > rlim, eff(R,r) = colleff(R,r) = 1). */

public:
  LongHydroProb()
      : hydroprob(SimmelTerminalVelocity{}), coaleff(1.0) {}

  LongHydroProb(const double coaleff)
      : hydroprob(SimmelTerminalVelocity{}), coaleff(coaleff) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns the probability of collision-coalescence
  using Simmel et al. 2002's formulation of Long's
  hydrodynamic, ie. gravitational, Kernel */
  {
    const auto eff = kerneleff(drop1, drop2);
    return hydroprob(drop1, drop2, eff, DELT, VOLUME);
  }
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_INLINE_FUNCTION
double LongHydroProb::kerneleff(const Superdrop &drop1,
                                const Superdrop &drop2) const
/* returns the efficiency of collision-coalescence, eff, according to
  equations 12 and 13 of Simmel et al. 2002). eff = eff(R,r) where R>r.
  eff = colleff(R,r) * coaleff(R,r). Usually it's assumed that
  coaleff(R,r) = 1, ie. eff = colleff, which also means that for
  collisions where R > rlim, eff(R,r) = colleff(R,r) = 1). */
{
  constexpr double rlim = 5e-5 / dlc::R0;          // 50 micron limit to determine collision-coalescence efficiency (eff)
  constexpr double colleff_lim = 0.001;            // minimum efficiency if larger droplet's radius < rlim
  constexpr double A1 = 4.5e4 * dlc::R0 * dlc::R0; // constants in efficiency calc if larger droplet's radius < rlim
  constexpr double A2 = 3e-4 / dlc::R0;

  const auto smallr = double{Kokkos::fmin(drop1.get_radius(), drop2.get_radius())};
  const auto bigr = double{Kokkos::fmax(drop1.get_radius(), drop2.get_radius())};

  /* calculate collision-coalescence efficiency, eff = colleff * coaleff */
  auto colleff = double{1.0};
  if (bigr < rlim)
  {
    const auto colleff_calc = double{A1 * bigr * bigr * (1 - A2 / smallr)};
    colleff = Kokkos::fmax(colleff_calc, colleff_lim); // colleff >= colleff_lim
  }

  const auto eff = colleff * coaleff;

  return eff;
}

#endif // LONGHYDROPROB_HPP