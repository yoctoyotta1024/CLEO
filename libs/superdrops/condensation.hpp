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
#include "impliciteuler.hpp"
#include "kokkosaliases_sd.hpp"
#include "microphysicalprocess.hpp"
#include "sdmmonitor.hpp"
#include "state.hpp"
#include "superdrop.hpp"
#include "thermodynamic_equations.hpp"

namespace dlc = dimless_constants;

/*
SuperdropletsChangeFunctor struct encapsulates superdroplet change during condensation
so that parallel loop in superdroplets_change function (see below) only captures
necessary objects and not other members of DoCondensation coincidentally
*/
struct SuperdropletsChangeFunctor {
  const ImplicitEuler &impe;    /**< Instance of ImplicitEuler ODE solver */
  const subviewd_supers supers; /** The view of superdroplets. */
  const State &state;           /**< The thermodynamic state. */
  const double s_ratio;         /**< s_ratio The saturation ratio. */
  const double ffactor;         /**< The sum of the diffusion factors. */

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

  /*
   * operator for functor in superdroplets_change function used in parallel (TeamThreadRangePolicy)
   * loop over superdroplets in supers view in order to call superdrop_mass_change
   */
  KOKKOS_INLINE_FUNCTION void operator()(const size_t kk, double &mass_condensed) const {
    const auto deltamass = superdrop_mass_change(supers(kk), state.temp, s_ratio, ffactor);
    mass_condensed += deltamass;
  }
};

/*
EffectOnThermodynamicStateFunctor struct encapsulates state change during condensation
so that parallel function effect_on_thermodynamic_state (see below) only captures
necessary objects and not other members of DoCondensation coincidentally
*/
struct EffectOnThermodynamicStateFunctor {
  const double totmass_condensed; /**< change in liquid mass in parcel volume 'dm' */

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

  /*
   * operator for functor in effect_on_thermodynamic_state function called in
   * parallel (Single PerTeam Policy) in order to call state_change
   */
  KOKKOS_INLINE_FUNCTION void operator()(State &state) const {
    constexpr double R0cubed_VOL0 = dlc::R0 * dlc::R0 * dlc::R0 / dlc::VOL0;
    const auto totrho_condensed = double{totmass_condensed / state.get_volume() *
                                         R0cubed_VOL0};  // drho_condensed_vapour/dt * delta t
    state = state_change(totrho_condensed, state);
  }
};

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
   * @param mo Monitor of SDM processes.
   */
  KOKKOS_FUNCTION
  void do_condensation(const TeamMember &team_member, const subviewd_supers supers, State &state,
                       const SDMMonitor auto mo) const {
    /* superdroplet radii changes */
    const auto totmass_condensed = superdroplets_change(team_member, supers, state);

    /* resultant effect on thermodynamic state */
    effect_on_thermodynamic_state(team_member, totmass_condensed, state);

    mo.monitor_condensation(team_member, totmass_condensed);
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
  KOKKOS_FUNCTION double superdroplets_change(const TeamMember &team_member,
                                              const subviewd_supers supers,
                                              const State &state) const;

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

 public:
  /**
   * @brief Constructs a DoCondensation object.
   * @param do_alter_thermo Whether to alter the thermodynamics of the State.
   * @param delt Time step to integrate ODE using implcit Euler method.
   * @param maxniters Maximum no. iterations of Newton Raphson Method.
   * @param rtol Relative tolerance for implicit Euler method.
   * @param atol Absolute tolerance for implicit Euler method.
   * @param minsubdelt Minimum subtimestep in cases of substepping implicit Euler method.
   */
  DoCondensation(const bool do_alter_thermo, const double delt, const size_t maxniters,
                 const double rtol, const double atol, const double minsubdelt)
      : do_alter_thermo(do_alter_thermo), impe(delt, maxniters, rtol, atol, minsubdelt) {}

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
   * @param mo Monitor of SDM processes.
   * @return The updated view super-droplets.
   */
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember &team_member, const unsigned int subt,
                                         subviewd_supers supers, State &state,
                                         const SDMMonitor auto mo) const {
    do_condensation(team_member, supers, state, mo);
  }
};

/**
 * @brief Constructs a microphysical process for condensation / evaporation of super-droplets
 * with a constant time-step 'interval'.
 *
 * @param interval The constant time-step for condensation.
 * @param do_alter_thermo Whether to alter the thermodynamic state after condensation / evaporation.
 * @param maxniters Maximum no. iterations of Newton Raphson Method.
 * @param step2dimlesstime A function to convert 'interval' time-step to a dimensionless time.
 * @param rtol Relative tolerance for implicit Euler method.
 * @param atol Absolute tolerance for implicit Euler method.
 * @param MINSUBDELT Minimum subtimestep in cases of substepping implicit Euler method.
 * @param realtime2dimless A function to convert a real-time to a dimensionless time.
 * @return The constructed microphysical process for condensation / evaporation.
 */
inline MicrophysicalProcess auto Condensation(
    const unsigned int interval, const std::function<double(unsigned int)> step2dimlesstime,
    const bool do_alter_thermo, const size_t maxniters, const double rtol, const double atol,
    const double MINSUBDELT, const std::function<double(double)> realtime2dimless) {
  const auto delt = step2dimlesstime(interval);  // dimensionless time equivalent to interval
  const auto minsubdelt = realtime2dimless(MINSUBDELT);  // dimensionless SUBDELT [s]

  const MicrophysicsFunc auto do_cond =
      DoCondensation(do_alter_thermo, delt, maxniters, rtol, atol, minsubdelt);

  return ConstTstepMicrophysics(interval, do_cond);
}

#endif  // LIBS_SUPERDROPS_CONDENSATION_HPP_
