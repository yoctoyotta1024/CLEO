/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: constprob.hpp
 * Project: collisions
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * constant probability of collisio event between two (real) droplets. Calculations is
 * contained in structure that satisfies the requirements of the PairProbability concept
 * (see collisions.hpp)
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_CONSTPROB_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_CONSTPROB_HPP_

#include <Kokkos_Core.hpp>

#include "../superdrop.hpp"

/* Probability of collision of a pair of droplets
as formulated in Shima et al. 2009 equation 3,
prob_jk = K(drop1, drop2) * delta_t/delta_vol.
Here  K(drop1, drop2) = K is a constant, e.g.
K = c + b where c and b are the rate of
coalescence and breakup (no rebound) */
struct ConstProb {
 private:
  double kernel;

 public:
  explicit ConstProb(const double k) : kernel(k) {}

  /* returns probability that a pair of droplets collide
  according to constant collision kernel, K(drop1, drop2) = const.
  Prob equation is : prob_jk = K(drop1, drop2) * delta_t/delta_vol */
  KOKKOS_INLINE_FUNCTION
  double operator()(const Superdrop &drop1, const Superdrop &drop2, const double DELT,
                    const double VOLUME) const {
    return kernel * DELT / VOLUME;  // prob_jk * time interval / volume [s/m^3]
  }
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_CONSTPROB_HPP_
