/*
 * ----- CLEO -----
 * File: golovinprob.hpp
 * Project: collisionprobs
 * Created Date: Thursday 9th November 2023
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
 * Header file for calculation of probability of a
 * collision-coalescence event between two (real)
 * droplets using the Golovin Kernel. Probability
 * calculations are contained in structures
 * that satisfy the requirements of the
 * PairProbability concept (see collisions.hpp)
*/

#ifndef GOLOVINPROB_HPP
#define GOLOVINPROB_HPP

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"

namespace dlc = dimless_constants;

struct GolovinProb
/* Probability of collision-coalescence of
a pair of droplets according to Golovin 1963
(see e.g. Shima et al. 2009) */
{
private:
  const double prob_jk_const;

public:
  GolovinProb()
      : prob_jk_const(1.5e3 * dlc::R0 * dlc::R0 * dlc::R0) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1,
                    const Superdrop &drop2,
                    const double DELT,
                    const double VOLUME) const;
  /* returns probability that a pair of droplets coalesces
  according to Golovin's (sum of volumes) coalescence kernel.
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where
  K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3),
  and K(drop1, drop2) is Golovin 1963 (coalescence) kernel */
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_INLINE_FUNCTION
double GolovinProb::operator()(const Superdrop &drop1,
                               const Superdrop &drop2,
                               const double DELT,
                               const double VOLUME) const
/* returns probability that a pair of droplets coalesces
according to Golovin's (sum of volumes) coalescence kernel.
Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where
K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3),
and K(drop1, drop2) is Golovin 1963 (coalescence) kernel */
{
  const auto DELT_DELVOL = double{DELT / VOLUME}; // time interval / volume for which collision probability is calculated [s/m^3]
  const auto golovins_kernel = double{prob_jk_const *
                                      (drop1.vol() + drop2.vol())}; // Golovin 1963 coalescence kernel

  const auto prob_jk = golovins_kernel * DELT_DELVOL;

  return prob_jk;
}

#endif // GOLOVINPROB_HPP
