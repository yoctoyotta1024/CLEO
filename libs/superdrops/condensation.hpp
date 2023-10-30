/*
 * ----- CLEO -----
 * File: condensation.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 30th October 2023
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
#include <math.h> 

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "./impliciteuler.hpp"
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
  bool doAlterThermo;          // whether to make condensation alter ThermoState or not
  ImplicitEuler impe;                 // implicit euler solver

  KOKKOS_FUNCTION
  void do_condensation(const subviewd_supers supers, State &state) const;
  /* Enacts condensation / evaporation microphysical process.
  Change to superdroplet radii and temp, qv and qc due to
  sum of radii changes via diffusion and condensation of
  water vapour during timestep delt. Using equations
  from "An Introduction To Clouds...." (see note at top of file) */

  KOKKOS_FUNCTION
  double condensation_mass_change(Superdrop &drop,
                                  const double temp,
                                  const double s_ratio,
                                  const double ffactor) const;

  KOKKOS_FUNCTION
  void condensation_state_change(const double totrho_condensed,
                                 State &state) const;
  /* change the thermodynamic variables (temp, qv and qc) of
  ThermoState state given the total change in condensed
  water mass per volume during time interval delt */

public:
  DoCondensation(const bool doAlterThermo,
                 const unsigned int niters,
                 const double delt,
                 const double maxrtol,
                 const double maxatol,
                 const double subdelt)
      : doAlterThermo(doAlterThermo),
        impe(niters, delt, maxrtol, maxatol, subdelt) {}

  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION
      subviewd_supers
      operator()(const unsigned int subt,
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
Condensation(const unsigned int interval,
             const bool doAlterThermo,
             const unsigned int niters,
             const std::function<double(int)> step2dimlesstime,
             const double maxrtol,
             const double maxatol,
             const double SUBDELT,
             const std::function<double(int)> realtime2dimless)

/* constructs Microphysical Process for
condensation/evaporation of superdroplets with a
constant timestep 'interval' given the
"do_condensation" function-like type */
{
  const double delt = step2dimlesstime(interval);   // dimensionless time [] equivlent to interval
  const double subdelt = realtime2dimless(SUBDELT); // dimensionless time [] equivlent to SUBDELT [s]

  const auto do_cond = DoCondensation(doAlterThermo, niters,
                                      delt, maxrtol, maxatol,
                                      subdelt);

  return ConstTstepMicrophysics(interval, do_cond);
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
  const double press(state.press);
  const double temp(state.temp);
  const double qvap(state.qvap);

  /* superdroplet radii changes */
  constexpr double C0cubed(dlc::COORD0 * dlc::COORD0 * dlc::COORD0);
  const double VOLUME(state.get_volume() * C0cubed);    // volume in which condensation occurs [m^3]
  const double psat(saturation_pressure(temp));
  const double s_ratio(supersaturation_ratio(press, qvap, psat));
  const double ffactor(diffusion_factor(press, temp, psat));

  double totmass_condensed(0.0); // cumulative change to liquid mass in parcel volume 'dm'
  for (size_t kk(0); kk < supers.extent(0); ++kk)
  {
    const double deltamass_condensed(
        condensation_mass_change(supers(kk), temp, s_ratio, ffactor));
    totmass_condensed += deltamass_condensed; // dm += dm_condensed_vapour/dt * delta t
  }
  const double totrho_condensed(totmass_condensed / VOLUME); // drho_condensed_vapour/dt * delta t

  /* resultant effect on thermodynamic state */
  if (doAlterThermo)
  {
    condensation_state_change(totrho_condensed, state);
  }
}

KOKKOS_FUNCTION
double DoCondensation::condensation_mass_change(Superdrop &drop,
                                                const double temp,
                                                const double s_ratio,
                                                const double ffactor) const
/* update superdroplet radius due to radial growth/shrink
  via condensation and diffusion of water vapour according
  to equations from "An Introduction To Clouds...." (see
  note at top of file). Then return mass of liquid that
  condensed onto /evaporated off of droplet. New radius is
  calculated using impliciteuler method which iterates
  condensation-diffusion ODE given the previous radius. */
{
  /* do not pass r by reference here!! copy value into iterator */
  const auto ab_kohler = kohler_factors(drop, temp); // pair = {akoh, bkoh}
  const double newr(impe.solve_condensation(s_ratio, ab_kohler, ffactor,
                                            drop.get_radius())); // timestepping eqn [7.28] forward
  const double delta_radius(drop.change_radius(newr));

  constexpr double R0cubed = dlc::R0 * dlc::R0 * dlc::R0;
  constexpr double dmdt_const = 4.0 * M_PI * dlc::Rho_l * R0cubed;
  const double rsqrd(drop.get_radius() * drop.get_radius());
  const double mass_condensed = (dmdt_const * rsqrd *
                                 drop.get_xi() * delta_radius); // eqn [7.22] * delta t

  return mass_condensed;
}

KOKKOS_FUNCTION void
DoCondensation::condensation_state_change(const double totrho_condensed,
                                         State &state) const
/* change the thermodynamic variables (temp, qv and qc) of
ThermoState state given the total change in condensed
water mass per volume during time interval delt */
{
  const double delta_qcond = totrho_condensed / dlc::Rho_dry;
  const double delta_temp = (dlc::Latent_v /
                             moist_specifc_heat(state.qvap, state.qcond)) *
                            delta_qcond;

  state.temp += delta_temp;
  state.qvap -= delta_qcond;
  state.qcond += delta_qcond;
}

#endif // CONDENSATION_HPP