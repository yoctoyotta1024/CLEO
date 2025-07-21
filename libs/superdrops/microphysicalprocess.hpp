/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: microphysicalprocess.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Tobias Kölling (TK)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Concept of a Microphysical Process as well as helpers structs and functions for creating
 * structs that obeys the concept to model microphysics in SDM (eg. condensation or
 * collision-coalescence using the ConstTstepMicrophysics struct).
 */

#ifndef LIBS_SUPERDROPS_MICROPHYSICALPROCESS_HPP_
#define LIBS_SUPERDROPS_MICROPHYSICALPROCESS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <concepts>

#include "../cleoconstants.hpp"
#include "kokkosaliases_sd.hpp"
#include "sdmmonitor.hpp"
#include "state.hpp"
#include "superdrop.hpp"

/**
 * @brief Concept of a microphysical process.
 *
 * The MicrophysicalProcess concept represents all types that meet the requirements (constraints)
 * of two time-stepping functions ("next_step" and "on_step"), as well as the constraints on the
 * "run_step" function.
 *
 * Note: NullSDMMonitor used here as placeholder for templated run_step function that can take any
 * type satisfying the SDMMonitor concept.
 *
 * @tparam P The type that satisfies the MicrophysicalProcess concept.
 */
template <typename P>
concept MicrophysicalProcess =
    requires(P p, const TeamMember &tm, const unsigned int t, subviewd_supers supers, State &state,
             const NullSDMMonitor mo) {
      { p.next_step(t) } -> std::convertible_to<unsigned int>;
      { p.on_step(t) } -> std::same_as<bool>;
      { p.run_step(tm, t, supers, state, mo) } -> std::same_as<void>;
    };

/**
 * @brief Combined microphysical process struct.
 *
 * The CombinedMicrophysicalProcess struct combines two microphysical processes into one.
 * It implements the MicrophysicalProcess concept by delegating calls to the individual processes.
 * Structure enacts associative addition operation that defines the set for the
 * microphysical process Monoid.
 *
 * @tparam Microphys1 The type of the first microphysical process.
 * @tparam Microphys2 The type of the second microphysical process.
 */
template <MicrophysicalProcess Microphys1, MicrophysicalProcess Microphys2>
struct CombinedMicrophysicalProcess {
 private:
  Microphys1 a; /**< The first instance of type of MicrophysicalProcess. */
  Microphys2 b; /**< The second instance of type of MicrophysicalProcess. */

 public:
  /**
   * @brief Constructs a CombinedMicrophysicalProcess object.
   *
   * @param a The first microphysical process.
   * @param b The second microphysical process.
   */
  CombinedMicrophysicalProcess(const Microphys1 a, const Microphys2 b) : a(a), b(b) {}

  /**
   * @brief Returns the next time step for the combined microphysical process.
   *
   * @param subt The current time step.
   * @return The smaller of the next time steps from the two individual processes.
   */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int subt) const {
    const auto t_a = a.next_step(subt);
    const auto t_b = b.next_step(subt);

    return Kokkos::min(t_a, t_b);
  }

  /**
   * @brief Checks if the combined microphysical process should perform an on-step action.
   *
   * @param subt The current time step.
   * @return True if either individual process indicates an on-step action.
   */
  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int subt) const { return a.on_step(subt) || b.on_step(subt); }

  /**
   * @brief Runs the combined microphysical process.
   *
   * @param team_member The Kokkos team member executing the process.
   * @param subt The current time step.
   * @param supers The view of super-droplets.
   * @param state The state of the system / volume.
   * @param mo Monitor of SDM processes.
   * @return The updated view of super-droplets after the process.
   */
  KOKKOS_INLINE_FUNCTION void run_step(const TeamMember &team_member, const unsigned int subt,
                                       subviewd_supers supers, State &state,
                                       const SDMMonitor auto mo) const {
    a.run_step(team_member, subt, supers, state, mo);
    b.run_step(team_member, subt, supers, state, mo);
  }
};

/**
 * @brief Operator for combining two microphysical processes.
 *
 * This operator combines two microphysical processes into one using the
 * CombinedMicrophysicalProcess struct.
 *
 * @param a The first microphysical process.
 * @param b The second microphysical process.
 * @return The combined microphysical process.
 */
auto operator>>(const MicrophysicalProcess auto a, const MicrophysicalProcess auto b) {
  return CombinedMicrophysicalProcess(a, b);
}

/**
 * @brief Null microphysical process struct.
 *
 * The NullMicrophysicalProcess struct represents a microphysical process that does nothing.
 * It is defined to satisfy null member of the Monoid set.
 */
struct NullMicrophysicalProcess {
  /**
   * @brief Returns the next time step for the null microphysical process.
   *
   * @param subt The current time step.
   * @return The maximum unsigned integer value, indicating no further time step will require
   * action of null microphyisical process.
   */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int subt) const { return LIMITVALUES::uintmax; }

  /**
   * @brief Checks if the null microphysical process should perform an on-step action.
   *
   * @param subt The current time step.
   * @return Always returns false, indicating no action is ever performed.
   */
  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int subt) const { return false; }

  /**
   * @brief Runs the null microphysical process.
   *
   * Is null, i.e. does nothing and returns unchanged super-droplet view.
   *
   * @param team_member The team member executing the process.
   * @param subt The current time step.
   * @param supers The view of super-droplets.
   * @param state The state of the system.
   * @param mo Monitor of SDM processes.
   * @return The unchanged view of super-droplets.
   */
  KOKKOS_INLINE_FUNCTION void run_step(const TeamMember &team_member, const unsigned int subt,
                                       subviewd_supers supers, State &state,
                                       const SDMMonitor auto mo) const {}
};

/**
 * @brief Concept for a microphysics function.
 *
 * The MicrophysicsFunc concept represents all function-like types that can be called by the
 * "run_step" function in ConstTstepMicrophysics.
 *
 * Note: NullSDMMonitor used here as placeholder for templated run_step function that can take any
 * type satisfying the SDMMonitor concept.
 *
 * @tparam F The type that satisfies the MicrophysicsFunc concept.
 */
template <typename F>
concept MicrophysicsFunc = requires(F f, const TeamMember &tm, const unsigned int subt,
                                    subviewd_supers supers, State &state, const NullSDMMonitor mo) {
  { f(tm, subt, supers, state, mo) } -> std::same_as<void>;
};

/**
 * @brief Struct representing microphysics with constant time step.
 *
 * The ConstTstepMicrophysics struct is a type that satisfies the concept of microphysical process
 * and has a constant time step interval. It can be used to create microphysical processes with
 * constant time steps between action of microphysics determined by the MicrophysicsFunc type 'F'.
 *
 * Special case: If interval is largest possible unsigned integer, on_step never returns true.
 *
 * @tparam F The type of the microphysics function.
 */
template <MicrophysicsFunc F>
struct ConstTstepMicrophysics {
 private:
  unsigned int interval; /**< The constant time step between calls to microphysics. */
  F do_microphysics;     /**< Function-like microphysics is type of MicrophysicsFunc*/

 public:
  /**
   * @brief Constructs a ConstTstepMicrophysics object.
   *
   * @param interval The constant time step between calls to microphysics.
   * @param f The microphysics function.
   */
  ConstTstepMicrophysics(const unsigned int interval, const F f)
      : interval(interval), do_microphysics(f) {}

  /**
   * @brief Returns the next time when the microphysics should be called given its constant
   * interval.
   *
   * @param subt The current time step.
   * @return The next time step based on the interval.
   */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int subt) const {
    return ((subt / interval) + 1) * interval;
  }

  /**
   * @brief Checks if the constant time step microphysics should perform an on-step action.
   *
   * Special case: If interval is largest possible unsigned integer, on_step never returns true.
   *
   * @param subt The current time step.
   * @return True if the current time step is a multiple of the interval.
   */
  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int subt) const {
    return (subt % interval == 0) && (interval != LIMITVALUES::uintmax);
  }

  /**
   * @brief Runs microphysics with the constant time step.
   *
   * @param team_member The team member executing the process.
   * @param subt The current time step.
   * @param supers The view of super-droplets.
   * @param state The state of the system / volume.
   * @param mo Monitor of SDM processes.
   * @return The updated view of super-droplets after the process.
   */
  KOKKOS_INLINE_FUNCTION void run_step(const TeamMember &team_member, const unsigned int subt,
                                       subviewd_supers supers, State &state,
                                       const SDMMonitor auto mo) const {
    if (on_step(subt)) {
      do_microphysics(team_member, subt, supers, state, mo);
    }
  }
};

#endif  // LIBS_SUPERDROPS_MICROPHYSICALPROCESS_HPP_
