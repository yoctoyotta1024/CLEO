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

struct LowListProb
/* Probability of collision-coalescence of a pair of
droplets as formulated in Shima et al. 2009 equation 3,
prob_jk = K(drop1, drop2) * delta_t/delta_vol.
Here K(drop1, drop2) is the hydrodynamic kernel 
with efficiency, eff = colleff * coaleff and
- colleff = Long's collision efficiency as
  given by equation 13 of Simmel et al. 2002.
- coaleff = Low and List 1982(a), equations (4.5) and (4.6)
*/
{
private:
  LongHydroProb longprob;

  KOKKOS_INLINE_FUNCTION
  double coaleff(const Superdrop &drop1,
                 const Superdrop &drop2) const;
  /* returns the coalescence efficiency according to
  Low and List 1982(a), equations (4.5) and (4.6)*/

public:
  LowListProb() : longprob() {}

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


#endif // LOWLISTPROB_HPP

