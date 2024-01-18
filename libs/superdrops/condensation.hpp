/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: condensation.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 17th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
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
#include <Kokkos_MathematicalConstants.hpp> // for pi
#include <Kokkos_Random.hpp>

#include "../cleoconstants.hpp"
#include "./impliciteuler.hpp"
#include "./kokkosaliases_sd.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"
#include "./state.hpp"
#include "./thermodynamic_equations.hpp"

namespace dlc = dimless_constants;

struct DoCondensation
/* function-like type that enacts
condensation / evaporation microphysical process */
{
private:
  bool doAlterThermo; // whether to make condensation alter ThermoState or not
  ImplicitEuler impe; // implicit euler solver

  KOKKOS_FUNCTION
  void do_condensation(const TeamMember &team_member,
                       const subviewd_supers supers,
                       State &state) const;
  /* Enacts condensation / evaporation microphysical process.
  Change to superdroplet radii and temp, qv and qc due to
  sum of radii changes via diffusion and condensation of
  water vapour during timestep delt. Using equations
  from "An Introduction To Clouds...." (see note at top of file) */

  KOKKOS_FUNCTION double
  superdroplets_change(const TeamMember &team_member,
                       const subviewd_supers supers,
                       const State &state) const;
  /* returns total change in liquid water mass
  in parcel volume 'mass_condensed' by enacting
  superdroplets' condensation / evaporation */

  KOKKOS_FUNCTION
  double superdrop_mass_change(Superdrop &drop,
                               const double temp,
                               const double s_ratio,
                               const double ffactor) const;
  /* update superdroplet radius due to radial growth/shrink
  via condensation and diffusion of water vapour according
  to equations from "An Introduction To Clouds...." (see
  note at top of file). Then return mass of liquid that
  condensed onto /evaporated off of droplet. New radius is
  calculated using impliciteuler method which iterates
  condensation-diffusion ODE given the previous radius. */

  KOKKOS_FUNCTION
  void effect_on_thermodynamic_state(
      const TeamMember &team_member,
      const double totmass_condensed,
      State &state) const;
  /* if doAlterThermo isn't false, use a single team
  member to change the state due to the effect
  of condensation / evaporation */

  KOKKOS_FUNCTION
  State state_change(const double totrho_condensed, State &state) const;
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

  KOKKOS_INLINE_FUNCTION subviewd_supers
  operator()(const TeamMember &team_member,
             const unsigned int subt,
             subviewd_supers supers,
             State &state,
             GenRandomPool genpool) const
  /* this operator is used as an "adaptor" for using
  condensation as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    do_condensation(team_member, supers, state);
    return supers;
  }
};

inline MicrophysicalProcess auto
Condensation(const unsigned int interval,
             const bool doAlterThermo,
             const unsigned int niters,
             const std::function<double(unsigned int)> step2dimlesstime,
             const double maxrtol,
             const double maxatol,
             const double SUBDELT,
             const std::function<double(double)> realtime2dimless)

/* constructs Microphysical Process for
condensation/evaporation of superdroplets with a
constant timestep 'interval' given the
"do_condensation" function-like type */
{
  const auto delt = step2dimlesstime(interval);   // dimensionless time [] equivlent to interval
  const auto subdelt = realtime2dimless(SUBDELT); // dimensionless time [] equivlent to SUBDELT [s]

  const auto do_cond = DoCondensation(doAlterThermo, niters,
                                      delt, maxrtol, maxatol,
                                      subdelt);

  return ConstTstepMicrophysics(interval, do_cond);
}

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION
void DoCondensation::
    do_condensation(const TeamMember &team_member,
                    const subviewd_supers supers,
                    State &state) const
/* Enacts condensation / evaporation microphysical process.
Change to superdroplet radii and temp, qv and qc due to
sum of radii changes via diffusion and condensation of
water vapour during timestep delt. Using equations
from "An Introduction To Clouds...." (see note at top of file) */
{
  /* superdroplet radii changes */
  double totmass_condensed(superdroplets_change(team_member,
                                                supers,
                                                state));

  /* resultant effect on thermodynamic state */
  effect_on_thermodynamic_state(team_member, totmass_condensed, state);
}

KOKKOS_FUNCTION
double DoCondensation::
    superdroplets_change(const TeamMember &team_member,
                         const subviewd_supers supers,
                         const State &state) const
/* returns total change in liquid water mass in parcel
volume, 'mass_condensed', by enacting superdroplets'
condensation / evaporation. Kokkos::parallel_reduce([...])
is equivalent to summing deltamass over for loop:
for (size_t kk(0); kk < nsupers; ++kk) {[...]}
when in serial*/
{
  const auto nsupers = (size_t)supers.extent(0);

  const auto psat = saturation_pressure(state.temp);
  const auto s_ratio = supersaturation_ratio(state.press, state.qvap, psat);
  const auto ffactor = diffusion_factor(state.press, state.temp, psat);

  auto totmass_condensed = double{0.0}; // cumulative change to liquid mass in parcel volume 'dm'
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team_member, nsupers),
      [&, this](const size_t kk, double &mass_condensed)
      {
        const auto deltamass = superdrop_mass_change(supers(kk), state.temp,
                                                     s_ratio, ffactor);
        mass_condensed += deltamass;
      },
      totmass_condensed);

  return totmass_condensed;
}

KOKKOS_FUNCTION
double DoCondensation::superdrop_mass_change(Superdrop &drop,
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
  const auto newr = impe.solve_condensation(s_ratio, ab_kohler, ffactor,
                                            drop.get_radius()); // timestepping eqn [7.28] forward
  const auto delta_radius = double{drop.change_radius(newr)};

  constexpr double R0cubed = dlc::R0 * dlc::R0 * dlc::R0;
  constexpr double dmdt_const = 4.0 * Kokkos::numbers::pi * dlc::Rho_l * R0cubed;
  const auto rsqrd = double{drop.get_radius() * drop.get_radius()};
  const auto mass_condensed = double{dmdt_const * rsqrd *
                                     drop.get_xi() * delta_radius}; // eqn [7.22] * delta t

  return mass_condensed;
}

KOKKOS_FUNCTION
void DoCondensation::effect_on_thermodynamic_state(
    const TeamMember &team_member,
    const double totmass_condensed,
    State &state) const
/* if doAlterThermo isn't false, use a single team
member to change the state due to the effect
of condensation / evaporation */
{
  Kokkos::single(
      Kokkos::PerTeam(team_member),
      [&, this](State &state)
      {
  if (doAlterThermo)
  {
    const auto VOLUME = double{state.get_volume() * dlc::VOL0};       // volume in which condensation occurs [m^3]
    const auto totrho_condensed = double{totmass_condensed / VOLUME}; // drho_condensed_vapour/dt * delta t
    state = state_change(totrho_condensed, state);
  } },
      state);
}

KOKKOS_FUNCTION State
DoCondensation::state_change(const double totrho_condensed,
                             State &state) const
/* change the thermodynamic variables (temp, qv and qc) of
ThermoState state given the total change in condensed
water mass per volume during time interval delt */
{
  const auto delta_qcond = double{totrho_condensed / dlc::Rho_dry};
  const auto delta_temp = double{(dlc::Latent_v /
                                  moist_specifc_heat(state.qvap, state.qcond)) *
                                 delta_qcond};

  state.temp += delta_temp;
  state.qvap -= delta_qcond;
  state.qcond += delta_qcond;

  return state;
}

#endif // CONDENSATION_HPP
