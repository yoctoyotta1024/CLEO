/*
 * ----- CLEO -----
 * File: condensation.hpp
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
 * struct for condensation / evaporation of water
 * causing diffusional growth / shrinking of
 * droplets in SDM. Equations referenced as (eqn [X.YY])
 * are from "An Introduction To Clouds From The 
 * Microscale to Climate" by Lohmann, Luond
 * and Mahrt, 1st edition.
 */

#ifndef CONDENSATION_HPP
#define CONDENSATION_HPP

#include <concepts>

#include <Kokkos_Core.hpp>

#include "./kokkosaliases_sd.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"
#include "./state.hpp"
#include "./thermodynamic_equations.hpp"
#include "./urbg.hpp"

struct DoCondensation
/* function-like type that enacts
condensation / evaporation microphysical process */
{
private:

  KOKKOS_FUNCTION
  void do_condensation(const subviewd_supers supers, State &state) const;
  /* Enacts condensation / evaporation microphysical process.
  Change to superdroplet radii and temp, qv and qc due to
  sum of radii changes via diffusion and condensation of
  water vapour during timestep delt. Using equations
  from "An Introduction To Clouds...." (see note at top of file) */

public:
  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION
  subviewd_supers operator()(const unsigned int subt,
                             subviewd_supers supers,
                             State &state,
                             URBG<DeviceType> urbg) const
  /* this operator is used as an "adaptor" for using
  condensation as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    do_condensation(supers, state);

    return supers;
  }
};

inline MicrophysicalProcess auto
Condensation(const unsigned int interval)
/* constructs Microphysical Process for
condensation/evaporation of superdroplets with a
constant timestep 'interval' given the
"do_condensation" function-like type */
{
  return ConstTstepMicrophysics(interval, DoCondensation{});
}

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION
void DoCondensation::do_condensation(const subviewd_supers supers,
                                     State &state) const
/* Enacts condensation / evaporation microphysical process.
Change to superdroplet radii and temp, qv and qc due to
sum of radii changes via diffusion and condensation of
water vapour during timestep delt. Using equations
from "An Introduction To Clouds...." (see note at top of file) */
{
  constexpr double C0cubed(dlc::COORD0 * dlc::COORD0 * dlc::COORD0);

  const double psat(saturation_pressure(state.temp));
  // const double s_ratio(supersaturation_ratio(state.press, state.qvap, psat));

  /* superdroplet radii changes */
  double tot_rho_condensed(0.0); // cumulative change to liquid mass in parcel volume
  for (size_t kk(0); kk < supers.extent(0); ++kk)
  {
    const double delta_mass_condensed(psat); // TODO
    // const double delta_mass_condensed = superdroplet_growth_by_condensation(state.press, state.temp,
    //                                                                         psat, s_ratio, delt,
    //                                                                         impliciteuler, SDinGBx.superdrop);
    const double VOLUME(state.get_volume() * C0cubed);    // volume in which condensation occurs [m^3]
    tot_rho_condensed += (delta_mass_condensed / VOLUME); // drho_condensed_vapour/dt * delta t
  }

  // /* resultant effect on thermodynamic state */
  // if (doAlterThermo)
  // {
  //   condensation_alters_thermostate(state, tot_rho_condensed); // TODO
  // }
}

#endif // CONDENSATION_HPP