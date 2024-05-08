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
 * Last Modified: Wednesday 8th May 2024
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

#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalConstants.hpp>  // for pi
#include <Kokkos_Random.hpp>
#include <concepts>

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
  bool do_alter_thermo; /**< Whether to make condensation/evaporation alter State or not */
  ImplicitEuler impe;   /**< instance of ImplicitEuler ODE solver */

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
   * if do_alter_thermo is true, use a single team member to change the thermodynamics of the
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
   * @param do_alter_thermo Whether to alter the thermodynamics of the State.
   * @param niters Number of iterations of implicit Euler method.
   * @param delt Time step to integrate ODE using implcit Euler method.
   * @param maxrtol Maximum relative tolerance for implicit Euler method.
   * @param maxatol Maximum absolute tolerance for implicit Euler method.
   * @param subdelt Sub-time step size in implicit Euler method.
   */
  DoCondensation(const bool do_alter_thermo, const unsigned int niters, const double delt,
                 const double maxrtol, const double maxatol, const double subdelt)
      : do_alter_thermo(do_alter_thermo), impe(niters, delt, maxrtol, maxatol, subdelt) {}

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
   * @return The updated view super-droplets.
   */
  KOKKOS_INLINE_FUNCTION subviewd_supers operator()(const TeamMember &team_member,
                                                    const unsigned int subt, subviewd_supers supers,
                                                    State &state) const {
    do_condensation(team_member, supers, state);
    return supers;
  }
};

/**
 * @brief Constructs a microphysical process for condensation / evaporation of super-droplets
 * with a constant time-step 'interval'.
 *
 * @param interval The constant time-step for condensation.
 * @param do_alter_thermo Whether to alter the thermodynamic state after condensation / evaporation.
 * @param niters Number of iterations of implicit Euler method.
 * @param step2dimlesstime A function to convert 'interval' time-step to a dimensionless time.
 * @param maxrtol Maximum relative tolerance for implicit Euler method.
 * @param maxatol Maximum absolute tolerance for implicit Euler method.
 * @param SUBDELT The sub-time step of implicit Euler method.
 * @param realtime2dimless A function to convert a real-time to a dimensionless time.
 * @return The constructed microphysical process for condensation / evaporation.
 */
inline MicrophysicalProcess auto Condensation(
    const unsigned int interval, const std::function<double(unsigned int)> step2dimlesstime,
    const bool do_alter_thermo, const unsigned int niters, const double maxrtol,
    const double maxatol, const double SUBDELT,
    const std::function<double(double)> realtime2dimless) {
  const auto delt = step2dimlesstime(interval);    // dimensionless time equivlent to interval
  const auto subdelt = realtime2dimless(SUBDELT);  // dimensionless time equivlent to SUBDELT [s]

  const auto do_cond = DoCondensation(do_alter_thermo, niters, delt, maxrtol, maxatol, subdelt);

  return ConstTstepMicrophysics(interval, do_cond);
}

#endif  // LIBS_SUPERDROPS_CONDENSATION_HPP_
