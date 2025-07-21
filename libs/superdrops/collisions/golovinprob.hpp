/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: golovinprob.hpp
 * Project: collisions
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for calculation of probability of a
 * collision-coalescence event between two (real)
 * droplets using the Golovin Kernel. Probability
 * calculations are contained in structures
 * that satisfy the requirements of the
 * PairProbability concept (see collisions.hpp)
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_GOLOVINPROB_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_GOLOVINPROB_HPP_

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"

namespace dlc = dimless_constants;

/* Probability of collision-coalescence of
a pair of droplets according to Golovin 1963
(see e.g. Shima et al. 2009) */
struct GolovinProb {
 private:
  double prob_jk_const;

 public:
  GolovinProb() : prob_jk_const(1.5e3 * dlc::R0 * dlc::R0 * dlc::R0) {}

  /* returns probability that a pair of droplets coalesces
  according to Golovin's (sum of volumes) coalescence kernel.
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol where
  K(drop1, drop2) := C(drop1, drop2) * |v1−v2|, (see Shima 2009 eqn 3),
  and K(drop1, drop2) is Golovin 1963 (coalescence) kernel */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop1, const Superdrop &drop2, const double DELT,
                    const double VOLUME) const;
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_GOLOVINPROB_HPP_
