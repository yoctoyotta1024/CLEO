/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: hydrodynamicprob.hpp
 * Project: collisions
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Probability of some kind of collision  event between two (real) droplets
 * using the hydrodynamic (i.e. gravitational) kernel
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_HYDRODYNAMICPROB_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_HYDRODYNAMICPROB_HPP_

#include <Kokkos_Core.hpp>

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"
#include "../terminalvelocity.hpp"

namespace dlc = dimless_constants;

template <VelocityFormula TerminalVelocity>
struct HydrodynamicProb {
 private:
  double prob_jk_const;
  TerminalVelocity terminalv;

 public:
  explicit HydrodynamicProb(TerminalVelocity tv)
      : prob_jk_const(Kokkos::numbers::pi * dlc::R0 * dlc::R0 * dlc::W0), terminalv(tv) {}

  /* returns probability that a pair of droplets collide
  (and coalesce or breakup etc.) according to the hydrodynamic,
  i.e. gravitational, collision kernel. Probability is given by
  prob_jk = K(drop1, drop2) * delta_t/delta_vol, (see Shima 2009 eqn 3)
  where the kernel, K(drop1, drop2) := eff * pi * (r1 + r2)^2 * |v1−v2|,
  given the efficiency factor eff = eff(drop1, drop2), for
  example as expressed in equation 11 of Simmel at al. 2002 for
  collision-coalescence */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop1, const Superdrop &drop2, const double eff,
                    const double DELT, const double VOLUME) const {
    /* time interval / volume for which
    probability is calculated [s/m^3] */
    const auto DELT_DELVOL = double{DELT / VOLUME};

    /* calculate Hydrodynamic Kernel*/
    const auto sumr = double{drop1.get_radius() + drop2.get_radius()};
    const auto sumrsqrd = double{sumr * sumr};
    const auto vdiff = double{Kokkos::abs(terminalv(drop1) - terminalv(drop2))};
    const auto hydro_kernel = double{prob_jk_const * eff * sumrsqrd * vdiff};

    /* calculate probability prob_jk analogous Shima 2009 eqn 3 */
    const auto prob_jk = hydro_kernel * DELT_DELVOL;

    return prob_jk;
  }
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_HYDRODYNAMICPROB_HPP_
