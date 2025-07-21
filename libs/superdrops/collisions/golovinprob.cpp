/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: golovinprob.cpp
 * Project: collisions
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to calculate the probability of a
 * collision-coalescence event between two (real)
 * droplets using the Golovin Kernel. Probability
 * calculations are contained in structures
 * that satisfy the requirements of the
 * PairProbability concept (see collisions.hpp)
 */

#include "./golovinprob.hpp"

/* returns probability that a pair of droplets coalesces
according to Golovin's (sum of volumes) coalescence kernel.
Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where
K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3),
and K(drop1, drop2) is Golovin 1963 (coalescence) kernel */
KOKKOS_FUNCTION
double GolovinProb::operator()(const Superdrop &drop1, const Superdrop &drop2, const double DELT,
                               const double VOLUME) const {
  const auto DELT_DELVOL = double{
      DELT /
      VOLUME};  // time interval / volume for which collision probability is calculated [s/m^3]
  const auto golovins_kernel =
      double{prob_jk_const * (drop1.vol() + drop2.vol())};  // Golovin 1963 coalescence kernel

  const auto prob_jk = golovins_kernel * DELT_DELVOL;

  return prob_jk;
}
