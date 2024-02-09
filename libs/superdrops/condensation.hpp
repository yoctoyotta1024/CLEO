/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: condensation.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Shin-ichiro Shima (SiS)
 * -----
 * Last Modified: Friday 9th February 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct and functions for makign microphyiscal process that enacts condensation / evaporation of
 * water via diffusional growth / shrinking of droplets in SDM. Equations referenced as
 * (eqn [X.YY]) are from "An Introduction To Clouds From The Microscale to Climate" by Lohmann,
 * Luond and Mahrt, 1st edition.
 */

#ifndef LIBS_SUPERDROPS_CONDENSATION_HPP_
#define LIBS_SUPERDROPS_CONDENSATION_HPP_

#include <concepts>

#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalConstants.hpp>  // for pi
#include <Kokkos_Random.hpp>

#include "../cleoconstants.hpp"
#include "./impliciteuler.hpp"
#include "./kokkosaliases_sd.hpp"
#include "./microphysicalprocess.hpp"
#include "./state.hpp"
#include "./superdrop.hpp"
#include "./thermodynamic_equations.hpp"

namespace dlc = dimless_constants;

/**
 * @struct DoCondensation
 * @brief Implements condensation and evaporation microphysics for super-droplets.
 */
struct DoCondensation {
 private:
  bool doAlterThermo;   /**< Whether to make condensation alter State or not */
  ImplicitEuler impe;   /**< Implicit Euler solver */

  /**
   * @brief Enacts condensation / evaporation microphysics.
   *
   * Enacts condensation / evaporation microphysics. Change to superdroplet radius, and
   * optionally thermodynamics of the State due to sum of water condensed via diffusion and
   * condensation / evporation of water vapour during a given timestep delt. Using equations
   * (eqn [X.YY]) from "An Introduction To Clouds From The Microscale to Climate" by Lohmann,
   * Luond and Mahrt, 1st edition.
   *
   * @param team_member The Kokkos team member.
   * @param supers The superdroplets.
   * @param state The state.
   */
  KOKKOS_FUNCTION
  void do_condensation(const TeamMember &team_member, const subviewd_supers supers,
                       State &state) const;

  /**
   * @brief Changes super-droplet radii according to condensation / evaporation and returns the
   * total change in liquid water mass in volume as a result.
   *
   * returns total change in liquid water mass (dimensionless) in volume, 'mass_condensed', due to
   * condensation onto / evaporation of super-droplets.
   *
   * The equivalent serial version of Kokkos::parallel_reduce([...]) is summing deltamass over loop:
   * @code
   * for (size_t ii(0); ii < ngbxs; ++ii)
   * {
   *  [...]
   *  totmass_condensed += deltamass;
   * }
   * @endcode
   *
   * @param team_member The Kokkos team member.
   * @param supers The superdroplets.
   * @param state The state.
   * @return The total change in liquid water mass.
   */
  KOKKOS_FUNCTION double superdroplets_change(const TeamMember &team_member,
                                              const subviewd_supers supers,
                                              const State &state) const;

  /**
   * @brief Updates the super-droplet radius and returns the mass of liquid condensed or evaporated.
   *
   *
   * Updates the super-droplet radius due to radial growth/shrink via condensation and diffusion of
   * water vapour according to equations from "An Introduction To Clouds From The Microscale to
   * Climate" by Lohmann, Luond and Mahrt, 1st edition. New radius is calculated using 'impe'
   * ImplicitEuler instance which iteratively solves forward integration of condensation-diffusion
   * ODE. Return mass of liquid that condensed onto / evaporated off of droplet.
   *
   * @param drop The super-droplet.
   * @param temp The ambient temperature.
   * @param s_ratio The saturation ratio.
   * @param ffactor The sum of the diffusion factors.
   * @return The mass of liquid condensed or evaporated.
  */
  KOKKOS_FUNCTION
  double superdrop_mass_change(Superdrop &drop, const double temp, const double s_ratio,
                               const double ffactor) const;

  /**
   * @brief Applies the effect of condensation / evaporation on the thermodynamics of the State.
   *
   * if doAlterThermo is true, use a single team member to change the thermodynamics of the
   * State due to the effect of condensation / evaporation.
   *
   * @param team_member The Kokkos team member.
   * @param totmass_condensed The total mass of liquid condensed.
   * @param state The State of the volume containing the super-droplets
   * (prior to condensation / evaporation).
   */
  KOKKOS_FUNCTION
  void effect_on_thermodynamic_state(const TeamMember &team_member, const double totmass_condensed,
                                     State &state) const;

  /**
   * @brief Changes the thermodynamic variables of the State.
   *
   * Changes the thermodynamic variables, temperature, vapour and liquid mass mixing ratios
   * (qvap and qcond respectively) of the State given the total change in condensed water mass
   * in its volume.
   *
   * @param totrho_condensed The total condensed water mass in volume of State.
   * @param state The State of the volume containing the super-droplets
   * (prior to condensation / evaporation).
   * @return The updated State.
   */
  KOKKOS_FUNCTION
  State state_change(const double totrho_condensed, State &state) const;

 public:
  /**
   * @brief Constructs a DoCondensation object.
   * @param doAlterThermo Whether to alter the thermodynamics of the State.
   * @param niters Number of iterations of implicit Euler method.
   * @param delt Time step to integrate ODE using implcit Euler method.
   * @param maxrtol Maximum relative tolerance for implicit Euler method.
   * @param maxatol Maximum absolute tolerance for implicit Euler method.
   * @param subdelt Sub-time step size in implicit Euler method.
   */
  DoCondensation(const bool doAlterThermo, const unsigned int niters, const double delt,
                 const double maxrtol, const double maxatol, const double subdelt)
      : doAlterThermo(doAlterThermo), impe(niters, delt, maxrtol, maxatol, subdelt) {}

  /**
   * @brief Operator used as an "adaptor" for using condensation as the function-like type
   * satisfying the MicrophysicsFunction concept.
   *
   * This operator is an "adaptor" for using condensation as the MicrophysicsFunction type in a
   * ConstTstepMicrophysics instance (*hint* which satsifies the MicrophysicalProcess concept).
   *
   * @param team_member The Kokkos team member.
   * @param subt The microphysics time step.
   * @param supers The view of super-droplets.
   * @param state The State.
   * @param genpool The Kokkos thread-safe random number generator pool.
   * @return The updated view super-droplets.
   */
  KOKKOS_INLINE_FUNCTION subviewd_supers operator()(const TeamMember &team_member,
                                                    const unsigned int subt, subviewd_supers supers,
                                                    State &state, GenRandomPool genpool) const {
    do_condensation(team_member, supers, state);
    return supers;
  }
};

/**
 * @brief Constructs a microphysical process for condensation / evaporation of super-droplets
 * with a constant time-step 'interval'.
 *
 * @param interval The constant time-step for condensation.
 * @param doAlterThermo Whether to alter the thermodynamic state after condensation / evaporation.
 * @param niters Number of iterations of implicit Euler method.
 * @param step2dimlesstime A function to convert 'interval' time-step to a dimensionless time.
 * @param maxrtol Maximum relative tolerance for implicit Euler method.
 * @param maxatol Maximum absolute tolerance for implicit Euler method.
 * @param SUBDELT The sub-time step of implicit Euler method.
 * @param realtime2dimless A function to convert a real-time to a dimensionless time.
 * @return The constructed microphysical process for condensation / evaporation.
 */
inline MicrophysicalProcess auto Condensation(
    const unsigned int interval, const bool doAlterThermo, const unsigned int niters,
    const std::function<double(unsigned int)> step2dimlesstime, const double maxrtol,
    const double maxatol, const double SUBDELT,
    const std::function<double(double)> realtime2dimless) {
  const auto delt = step2dimlesstime(interval);    // dimensionless time equivlent to interval
  const auto subdelt = realtime2dimless(SUBDELT);  // dimensionless time equivlent to SUBDELT [s]

  const auto do_cond = DoCondensation(doAlterThermo, niters, delt, maxrtol, maxatol, subdelt);

  return ConstTstepMicrophysics(interval, do_cond);
}

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

/**
 * @brief Enacts condensation / evaporation microphysics.
 *
 * Enacts condensation / evaporation microphysics. Change to superdroplet radius, and
 * optionally thermodynamics of the State due to sum of water condensed via diffusion and
 * condensation / evporation of water vapour during a given timestep delt. Using equations
 * (eqn [X.YY]) from "An Introduction To Clouds From The Microscale to Climate" by Lohmann,
 * Luond and Mahrt, 1st edition.
 *
 * @param team_member The Kokkos team member.
 * @param supers The superdroplets.
 * @param state The state.
 */
KOKKOS_FUNCTION
void DoCondensation::do_condensation(const TeamMember &team_member, const subviewd_supers supers,
                                     State &state) const {
  /* superdroplet radii changes */
  double totmass_condensed(superdroplets_change(team_member, supers, state));

  /* resultant effect on thermodynamic state */
  effect_on_thermodynamic_state(team_member, totmass_condensed, state);
}

/**
 * @brief Changes super-droplet radii according to condensation / evaporation and returns the
 * total change in liquid water mass in volume as a result.
 *
 * returns total change in liquid water mass (dimensionless) in volume, 'mass_condensed', due to
 * condensation onto / evaporation of super-droplets.
 *
 * The equivalent serial version of Kokkos::parallel_reduce([...]) is summing deltamass over loop:
 * @code
 * for (size_t ii(0); ii < ngbxs; ++ii)
 * {
 *  [...]
 *  totmass_condensed += deltamass;
 * }
 * @endcode
 *
 * @param team_member The Kokkos team member.
 * @param supers The superdroplets.
 * @param state The state.
 * @return The total change in liquid water mass.
 */
KOKKOS_FUNCTION
double DoCondensation::superdroplets_change(const TeamMember &team_member,
                                            const subviewd_supers supers,
                                            const State &state) const {
  const auto nsupers = static_cast<size_t>(supers.extent(0));

  const auto psat = saturation_pressure(state.temp);
  const auto s_ratio = supersaturation_ratio(state.press, state.qvap, psat);
  const auto ffactor = diffusion_factor(state.press, state.temp, psat);

  auto totmass_condensed = double{0.0};  // cumulative change to liquid mass in parcel volume 'dm'
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team_member, nsupers),
      [&, this](const size_t kk, double &mass_condensed) {
        const auto deltamass = superdrop_mass_change(supers(kk), state.temp, s_ratio, ffactor);
        mass_condensed += deltamass;
      },
      totmass_condensed);

  return totmass_condensed;
}

/**
 * @brief Updates the super-droplet radius and returns the mass of liquid condensed or evaporated.
 *
 *
 * Updates the super-droplet radius due to radial growth/shrink via condensation and diffusion of
 * water vapour according to equations from "An Introduction To Clouds From The Microscale to
 * Climate" by Lohmann, Luond and Mahrt, 1st edition. New radius is calculated using 'impe'
 * ImplicitEuler instance which iteratively solves forward integration of condensation-diffusion
 * ODE. Return mass of liquid that condensed onto / evaporated off of droplet.
 *
 * @param drop The super-droplet.
 * @param temp The ambient temperature.
 * @param s_ratio The saturation ratio.
 * @param ffactor The sum of the diffusion factors.
 * @return The mass of liquid condensed or evaporated.
*/
KOKKOS_FUNCTION
double DoCondensation::superdrop_mass_change(Superdrop &drop, const double temp,
                                             const double s_ratio, const double ffactor) const {
  /* do not pass r by reference here!! copy value into iterator */
  const auto ab_kohler = kohler_factors(drop, temp);  // pair = {akoh, bkoh}
  const auto newr = impe.solve_condensation(s_ratio, ab_kohler, ffactor,
                                            drop.get_radius());  // timestepping eqn [7.28] forward
  const auto delta_radius = double{drop.change_radius(newr)};

  constexpr double R0cubed = dlc::R0 * dlc::R0 * dlc::R0;
  constexpr double dmdt_const = 4.0 * Kokkos::numbers::pi * dlc::Rho_l * R0cubed;
  const auto rsqrd = double{drop.get_radius() * drop.get_radius()};
  const auto mass_condensed =
      double{dmdt_const * rsqrd * drop.get_xi() * delta_radius};  // eqn [7.22] * delta t

  return mass_condensed;
}

/**
 * @brief Applies the effect of condensation / evaporation on the thermodynamics of the State.
 *
 * if doAlterThermo is true, use a single team member to change the thermodynamics of the
 * State due to the effect of condensation / evaporation.
 *
 * @param team_member The Kokkos team member.
 * @param totmass_condensed The total mass of liquid condensed.
 * @param state The State of the volume containing the super-droplets
 * (prior to condensation / evaporation).
 */
KOKKOS_FUNCTION
void DoCondensation::effect_on_thermodynamic_state(const TeamMember &team_member,
                                                   const double totmass_condensed,
                                                   State &state) const {
  Kokkos::single(
      Kokkos::PerTeam(team_member),
      [&, this](State &state) {
        if (doAlterThermo) {
          const auto VOLUME =
              double{state.get_volume() * dlc::VOL0};  // volume in which condensation occurs [m^3]
          const auto totrho_condensed =
              double{totmass_condensed / VOLUME};  // drho_condensed_vapour/dt * delta t
          state = state_change(totrho_condensed, state);
        }
      },
      state);
}

/**
 * @brief Changes the thermodynamic variables of the State.
 *
 * Changes the thermodynamic variables, temperature, vapour and liquid mass mixing ratios
 * (qvap and qcond respectively) of the State given the total change in condensed water mass
 * in its volume.
 *
 * @param totrho_condensed The total condensed water mass in volume of State.
 * @param state The State of the volume containing the super-droplets
 * (prior to condensation / evaporation).
 * @return The updated State.
 */
KOKKOS_FUNCTION State DoCondensation::state_change(const double totrho_condensed,
                                                   State &state) const {
  const auto delta_qcond = double{totrho_condensed / dlc::Rho_dry};
  const auto delta_temp =
      double{(dlc::Latent_v / moist_specifc_heat(state.qvap, state.qcond)) * delta_qcond};

  state.temp += delta_temp;
  state.qvap -= delta_qcond;
  state.qcond += delta_qcond;

  return state;
}

#endif  // LIBS_SUPERDROPS_CONDENSATION_HPP_
