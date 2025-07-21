/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrop_attrs.cpp
 * Project: superdrops
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality file for definition of a superdropet.
 * Equations referenced as (eqn [X.YY]) are from
 * "An Introduction To Clouds From The Microscale
 * to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 */

#include "superdrop_attrs.hpp"

/**
 * @brief Change the radius of droplet.
 *
 * Update droplet radius to larger out of new radius 'newr' or dry radius and return the
 * resultant change in radius = new radius - old radius. Prevents drops shrinking further once
 * they are size of dry radius.
 *
 * @param newr The new radius to set.
 * @return The change in radius.
 */
KOKKOS_FUNCTION
double SuperdropAttrs::change_radius(const double newr) {
  const auto oldradius = radius;

  /* if droplets are dry, do not shrink further */
  const auto dryr = dryradius();
  radius = Kokkos::fmax(newr, dryr);  // Kokkos equivalent to std::max() for floats (gpu compatible)

  /* return change in radius due to growth/shrinking of droplet */
  return radius - oldradius;
}

/**
 * @brief Get the total droplet mass.
 *
 * Calculates and returns total droplet mass = water + dry areosol.
 *
 * @return The total droplet mass.
 */
KOKKOS_FUNCTION
double SuperdropAttrs::mass() const {
  constexpr double massconst(4.0 / 3.0 * Kokkos::numbers::pi * dlc::Rho_l);  // 4/3 * pi * density
  const auto density_factor = double{1.0 - dlc::Rho_l / solute.rho_sol()};   // to account for msol

  auto mass = double{msol * density_factor};  // mass contribution of solute
  mass += massconst * rcubed();

  return mass;
}
