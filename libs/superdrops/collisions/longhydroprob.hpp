/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: longhydroprob.hpp
 * Project: collisions
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Probability of collision-coalescence event between two (real) droplets using the
 * hydrodynamic (i.e. gravitational) kernel according to Simmel et al. 2002's formulation
 * of Long's Hydrodynamic Kernel. Probability calculations are contained in structures
 * that satisfy the requirements of the PairProbability concept (see collisions.hpp)
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_LONGHYDROPROB_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_LONGHYDROPROB_HPP_

#include <Kokkos_Core.hpp>

#include "../../cleoconstants.hpp"
#include "../superdrop.hpp"
#include "../terminalvelocity.hpp"
#include "./hydrodynamicprob.hpp"

namespace dlc = dimless_constants;

/* Opertator returns the collision-coalescence probability
  given the efficiency factor, eff = eff(drop1, drop2),
  from Long's hydrodynamic kernel according to Simmel et al. 2002.
*/
struct LongHydroProb {
 private:
  HydrodynamicProb<SimmelTerminalVelocity> hydroprob;
  double coaleff;

  /* returns the efficiency of collision-coalescence, eff, according to
  equations 12 and 13 of Simmel et al. 2002). eff = eff(R,r) where R>r.
  eff = colleff(R,r) * coaleff(R,r). Usually it's assumed that
  coaleff(R,r) = 1, ie. eff = colleff, which also means that for
  collisions where R > rlim, eff(R,r) = colleff(R,r) = 1). */
  KOKKOS_FUNCTION
  double kerneleff(const Superdrop &drop1, const Superdrop &drop2) const;

 public:
  LongHydroProb() : hydroprob(SimmelTerminalVelocity{}), coaleff(1.0) {}

  explicit LongHydroProb(const double coaleff)
      : hydroprob(SimmelTerminalVelocity{}), coaleff(coaleff) {}

  /* returns the probability of collision-coalescence
  using Simmel et al. 2002's formulation of Long's
  hydrodynamic, ie. gravitational, Kernel */
  KOKKOS_FUNCTION
  double operator()(const Superdrop &drop1, const Superdrop &drop2, const double DELT,
                    const double VOLUME) const {
    const auto eff = kerneleff(drop1, drop2);
    return hydroprob(drop1, drop2, eff, DELT, VOLUME);
  }
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_LONGHYDROPROB_HPP_
