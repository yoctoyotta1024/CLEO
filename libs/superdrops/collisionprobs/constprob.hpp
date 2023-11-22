/*
 * ----- CLEO -----
 * File: constprob.hpp
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
 * constant probability of collisio event
 * between two (real) droplets. Calculations is
 * contained in structure that satisfies the
 * requirements of the PairProbability concept
 * (see collisions.hpp)
 */

#ifndef CONSTPROB_HPP
#define CONSTPROB_HPP

#include <Kokkos_Core.hpp>

#include "../superdrop.hpp"

struct ConstProb
/* Probability of collision of a pair of droplets
as formulated in Shima et al. 2009 equation 3,
prob_jk = K(drop1, drop2) * delta_t/delta_vol.
Here  K(drop1, drop2) = K is a constant, e.g.
K = c + b where c and b are the rate of
coalescence and breakup (no rebound) */
{
private:
  const double kernel;

public:
  ConstProb(const double k) : kernel(k) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const
  /* returns probability that a pair of droplets collide
  according to constant collision kernel, K(drop1, drop2) = const.
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol */
  {
    return kernel * DELT / VOLUME; // prob_jk * time interval / volume [s/m^3]
  }
};


#endif // CONSTPROB_HPP

