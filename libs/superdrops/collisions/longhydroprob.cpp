/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: longhydroprob.cpp
 * Project: collisions
 * Created Date: Wednesday 22nd November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for probability of collision-coalescence event between two (real) droplets using
 * the hydrodynamic (i.e. gravitational) kernel according to Simmel et al. 2002's formulation
 * of Long's Hydrodynamic Kernel. Probability calculations are contained in structures
 * that satisfy the requirements of the PairProbability concept (see collisions.hpp)
 */

#include "./longhydroprob.hpp"

/* returns the efficiency of collision-coalescence, eff, according to
  equations 12 and 13 of Simmel et al. 2002). eff = eff(R,r) where R>r.
  eff = colleff(R,r) * coaleff(R,r). Usually it's assumed that
  coaleff(R,r) = 1, ie. eff = colleff, which also means that for
  collisions where R > rlim, eff(R,r) = colleff(R,r) = 1). */
KOKKOS_FUNCTION
double LongHydroProb::kerneleff(const Superdrop &drop1, const Superdrop &drop2) const {
  constexpr double k1 = 4.5e4 * dlc::R0 * dlc::R0 * 100 * 100;
  constexpr double k2 = 3e-4 / dlc::R0 / 100;
  constexpr double rlim = 5e-5 / dlc::R0;
  constexpr double colleff_min = 0.001;

  const auto rsmall = double{Kokkos::fmin(drop1.get_radius(), drop2.get_radius())};
  const auto rbig = double{Kokkos::fmax(drop1.get_radius(), drop2.get_radius())};

  auto colleff = double{1.0};
  if (rbig < rlim) {
    colleff = k1 * rbig * rbig * (1 - (k2 / rsmall));
    colleff = Kokkos::fmax(colleff, colleff_min);
  }

  const auto eff = colleff * coaleff;

  return eff;
}
