/*
 * ----- CLEO -----
 * File: lowlistprob.hpp
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
 * probability of collision event
 * between two (real) droplets. Calculations is
 * contained in structure that satisfies the
 * requirements of the PairProbability concept
 * (see collisions.hpp)
 */

#ifndef LOWLISTPROB_HPP
#define LOWLISTPROB_HPP

#include <Kokkos_Core.hpp>

#include "../superdrop.hpp"
#include "./longhydroprob.hpp"
#include "./collisionkinetics.hpp"

struct LowListCoalProb
/* Probability of collision-coalescence of a pair of
droplets as formulated in Shima et al. 2009 equation 3,
prob_jk = K(drop1, drop2) * delta_t/delta_vol.
Here K(drop1, drop2) is the hydrodynamic kernel
with efficiency, eff = colleff * coaleff and
- colleff = Long's collision efficiency as
  given by equation 13 of Simmel et al. 2002.
- coaleff = Low and List 1982(a), equations (4.5)
  and (4.6) (see also McFarquhar 2004).
*/
{
private:
  LongHydroProb longprob;

  KOKKOS_INLINE_FUNCTION
  double coaleff(const Superdrop &drop1,
                 const Superdrop &drop2) const;
  /* returns coaleff, the coalescence efficiency
  of two droplets (given that they have collided)
  according to equations (4.5) and (4.6)
  Low and List 1982(a) */

  KOKKOS_INLINE_FUNCTION
  double expon(const double etot,
               const double surf_c) const;
  /* returns the exponential factor in eqn 4.5
  Low and List 1982(a) given the total collision energy,
  etot [J] and equivalent surface energy, surf_c [J] */

  KOKKOS_INLINE_FUNCTION
  double sizeratio_factor(const double r1,
                          const double r2) const;
  /* returns factor that takes into account the size
  ratio of droplets in eqn 4.5 Low and List 1982(a). */

public:
  LowListCoalProb() : longprob() {}

  KOKKOS_INLINE_FUNCTION
  double get_coaleff(const Superdrop &drop1,
                     const Superdrop &drop2) const
  {
    return coaleff(drop1, drop2);
  }

  KOKKOS_INLINE_FUNCTION
  double get_longprob(const Superdrop &drop1,
                      const Superdrop &drop2,
                      const double DELT,
                      const double VOLUME) const
  {
    return longprob(drop1, drop2, DELT, VOLUME);
  }

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns probability of collision-coalescence for a
  pair of droplets according to Long's formulation of the
  hydrodyanmic kernel for the collision probability modified
  by the coalescence efficiency from Low and List 1982(a). */
  {
    return longprob(drop1, drop2, DELT, VOLUME) * coaleff(drop1, drop2);
  }
};

struct LowListBuProb
/* Probability of collision-breakup of a pair of
droplets as formulated in Shima et al. 2009 equation 3,
prob_jk = K(drop1, drop2) * delta_t/delta_vol.
Here K(drop1, drop2) is the hydrodynamic kernel
with efficiency, eff = colleff * coaleff and
- colleff = Long's collision efficiency as
  given by equation 13 of Simmel et al. 2002.
- bueff = breakup efficiency, bueff = 1-coaleff,
  where coaleff is from equations (4.5) and
  (4.6) of  Low and List 1982(a)
  (see also McFarquhar 2004).
*/
{
private:
  LowListCoalProb ll;

public:
  LowListBuProb() : ll() {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns probability of collision-coalescence for a
  pair of droplets according to Long's formulation of the
  hydrodyanmic kernel for the collision probability modified
  by the coalescence efficiency from Low and List 1982(a). */
  {
    const double bueff(1.0 - ll.get_coaleff(drop1, drop2));
    const double longprob(ll.get_longprob(drop1, drop2,
                                          DELT, VOLUME));
    return longprob * bueff;
  }
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */
KOKKOS_INLINE_FUNCTION
double coaleff(const Superdrop &drop1,
               const Superdrop &drop2) const
/* returns coaleff, the coalescence efficiency
of two droplets (given that they have collided)
according to equations (4.5) and (4.6)
Low and List 1982(a) */
{
  constexpr double aconst = 0.778;
  constexpr double elim = 5e-6; // total energy limit [J]

  const double r1(drop1.get_radius());
  const double r2(drop2.get_radius());

  const double cke(collision_kinetic_energy(r1, r2,
                                            terminalv(drop1),
                                            terminalv(drop2)));
  const double surf_t(total_surfenergy(r1, r2)); // [J] S_t
  const double surf_c(coal_surfenergy(r1, r2));  // [J] S_c
  const double etot(cke + surf_t - surf_c);            // [J] total energy

  if (etot < elim)
  {
    const double radiiratio = sizeratio_factor(r1, r2);
    const double coaleff = aconst * radiiratio * expon(etot, surf_c);

    return coaleff;
  }
  else // coaleff = 0.0
  {
    return 0.0;
  }
}

KOKKOS_INLINE_FUNCTION
double expon(const double etot,
            const double surf_c) const
/* returns the exponential factor in eqn 4.5
Low and List 1982(a) given the total collision energy,
etot [J] and equivalent surface energy, surf_c [J] */
{
  constexpr double bconst = -2.62e6; // [J^-2]
  constexpr double sigma = 7.28e-2;  // [J/m^-2]

  const double exponent(bconst * sigma * etot * etot / surf_c);

  return Kokkos::exp(exponent);
}

KOKKOS_INLINE_FUNCTION
double sizeratio_factor(const double r1,
                        const double r2) const
/* returns factor that takes into account the size
ratio of droplets in eqn 4.5 Low and List 1982(a). */
{
  const double rsmall(Kokkos::fmin(r1, r2));
  const double rbig(Kokkos::fmax(r1, r2));
  const double alpha(1 + rsmall / rbig); // alpha = 1 + Ds/Dl

  return 1.0 / (alpha * alpha); // alpha^(-2)
}

#endif // LOWLISTPROB_HPP
