/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: lowlistprob.hpp
 * Project: collisions
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * probability of collision event between two (real) droplets accordin to Low and List.
 * Calculations is contained in structure that satisfies the requirements of the PairProbability
 * concept (see collisions.hpp)
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_LOWLISTPROB_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_LOWLISTPROB_HPP_

#include <Kokkos_Core.hpp>

#include "../superdrop.hpp"
#include "./collisionkinetics.hpp"
#include "./longhydroprob.hpp"

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
struct LowListCoalProb {
 private:
  LongHydroProb longprob;

  /* returns coaleff, the coalescence efficiency
  of two droplets (given that they have collided)
  according to equations (4.5) and (4.6)
  Low and List 1982(a) */
  KOKKOS_FUNCTION
  double coaleff(const Superdrop &drop1, const Superdrop &drop2) const;

  /* returns the exponential factor in eqn 4.5
  Low and List 1982(a) given the total collision energy,
  etot [J] and equivalent surface energy, surf_c [J] */
  KOKKOS_FUNCTION
  double expon(const double etot, const double surf_c) const;

  /* returns factor that takes into account the size
  ratio of droplets in eqn 4.5 Low and List 1982(a). */
  KOKKOS_FUNCTION
  double sizeratio_factor(const double r1, const double r2) const;

 public:
  LowListCoalProb() : longprob() {}

  KOKKOS_FUNCTION
  double get_coaleff(const Superdrop &drop1, const Superdrop &drop2) const {
    return coaleff(drop1, drop2);
  }

  KOKKOS_FUNCTION
  double get_longprob(const Superdrop &drop1, const Superdrop &drop2, const double DELT,
                      const double VOLUME) const {
    return longprob(drop1, drop2, DELT, VOLUME);
  }

  /* returns probability of collision-coalescence for a
  pair of droplets according to Long's formulation of the
  hydrodynamic kernel for the collision probability modified
  by the coalescence efficiency from Low and List 1982(a). */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop1, const Superdrop &drop2, const double DELT,
                    const double VOLUME) const {
    return longprob(drop1, drop2, DELT, VOLUME) * coaleff(drop1, drop2);
  }
};

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
struct LowListBuProb {
 private:
  LowListCoalProb ll;

 public:
  LowListBuProb() : ll() {}

  /* returns probability of collision-coalescence for a
  pair of droplets according to Long's formulation of the
  hydrodynamic kernel for the collision probability modified
  by the coalescence efficiency from Low and List 1982(a). */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop1, const Superdrop &drop2, const double DELT,
                    const double VOLUME) const {
    const auto bueff = double{1.0 - ll.get_coaleff(drop1, drop2)};
    const auto longprob = double{ll.get_longprob(drop1, drop2, DELT, VOLUME)};
    return longprob * bueff;
  }
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_LOWLISTPROB_HPP_
