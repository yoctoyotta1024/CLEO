/*
 * ----- CLEO -----
 * File: superdrop_attrs.cpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 26th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functionality file for definition of a superdropet.
 * Equations referenced as (eqn [X.YY]) are from
 * "An Introduction To Clouds From The Microscale
 * to Climate" by Lohmann, Luond and Mahrt, 1st edition.
 */

#include "./superdrop_attrs.hpp"

KOKKOS_FUNCTION double SuperdropAttrs::mass() const
/* returns total droplet mass = water + dry areosol  */
{
  constexpr double massconst(4.0 / 3.0 * M_PI * dlc::Rho_l); // 4/3 * pi * density
  const double density_factor(1.0 - dlc::Rho_l / solute.rho_sol()); // to account for msol

  const double rcubed(radius * radius * radius); // radius cubed
  double mass(msol * density_factor); // mass contribution of solute
  mass += massconst * rcubed;

  return mass;
}