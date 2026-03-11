/*
 * Copyright (c) 2026 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: coalescence_efficiency.cpp
 * Project: collisions
 * Created Date: Monday 9th March 2026
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * various parameterisations of the coalescence efficiency
 */

#include "./coalescence_efficiency.hpp"

/*
 * coalescence efficency given a collision occurs according to parameterisation from
 * Straub et al. 2010 section 3, equation 5
 * and Schlottke et al. 2010 section 4a equation 11
 * cke == collision_kinetic_energy
 * */
KOKKOS_FUNCTION
double coalescence_efficiency_straub2010(const Superdrop& drop1, const Superdrop& drop2,
                                         const double cke) {
  constexpr double beta = -1.15;

  const auto surf_c = coal_surfenergy(drop1.get_radius(),
                                      drop2.get_radius());  // [J] S_c
  const auto weber = double{cke / surf_c};
  const auto ecoal = double{Kokkos::exp(beta * weber)};

  return ecoal;
}
