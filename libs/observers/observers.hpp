/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: observers.hpp
 * Project: observers
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 21st May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer Concept and related structures for various ways
 * of observing (outputing data from) CLEO.
 * An example of an observer is one that outputs some data
 * from a Gridbox's State to a computer screen.
 */

#ifndef LIBS_OBSERVERS_OBSERVERS_HPP_
#define LIBS_OBSERVERS_OBSERVERS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./sdmmonitor.hpp"

/**
 * @brief Concept Observer is all types that have functions for timestepping and
 * observations functions as constrained here.
 *
 * @tparam Obs The type that satisfies the Observer concept.
 */
template <typename Obs>
concept Observer = requires(Obs obs, unsigned int t, const viewd_constgbx d_gbxs,
                            const viewd_constsupers totsupers) {
  { obs.next_obs(t) } -> std::convertible_to<unsigned int>;
  { obs.on_step(t) } -> std::same_as<bool>;
  { obs.before_timestepping(d_gbxs) } -> std::same_as<void>;
  { obs.after_timestepping() } -> std::same_as<void>;
  { obs.at_start_step(t, d_gbxs, totsupers) } -> std::same_as<void>;
  { obs.get_monitor_of_sdm_processes() } -> std::same_as<NullSDMMonitor>;
  // TODO(CB): return type for SDMMonitor concept
};

/**
 * @brief Structure CombinedObserver represents a new observer formed from combination of two
 * Observers 'a' and 'b'.
 *
 * @tparam Obs1 Type satisfying the Observer concept.
 * @tparam Obs2 Type satisfying the Observer concept.
 */
template <Observer Obs1, Observer Obs2>
struct CombinedObserver {
 private:
  Obs1 a; /**< First Observer. */
  Obs2 b; /**< Second Observer. */

 public:
  /**
   * @brief Construct a new CombinedObserver object.
   *
   * @param obs1 First Observer.
   * @param obs2 Second Observer.
   */
  CombinedObserver(const Obs1 obs1, const Obs2 obs2) : a(obs1), b(obs2) {}

  /**
   * @brief Run before timestepping for combination of 2 observers.
   *
   * Each observer is run sequentially.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    a.before_timestepping(d_gbxs);
    b.before_timestepping(d_gbxs);
  }

  /**
   * @brief Run after timestepping for combination of 2 observers.
   *
   * Each observer is run sequentially.
   */
  void after_timestepping() const {
    a.after_timestepping();
    b.after_timestepping();
  }

  /**
   * @brief Determine the next observation time for combination of 2 observers.
   *
   * For combination of 2 observers, the next observation time is the smaller out of the two
   * possible.
   *
   * @param t_mdl The unsigned int parameter.
   * @return unsigned int The next observation time.
   */
  unsigned int next_obs(const unsigned int t_mdl) const {
    const auto t_a = a.next_obs(t_mdl);
    const auto t_b = b.next_obs(t_mdl);

    /* return smaller of two unsigned ints (see std::min) */
    return !(t_a < t_b) ? t_b : t_a;
  }

  /**
   * @brief Check if on_step = true for combination of 2 observers.
   *
   * For combination of 2 observers, return on_step = true if either observer is on_step. Else
   * return false.
   *
   * @param t_mdl The unsigned int parameter.
   * @return bool True if on step, false otherwise.
   */
  bool on_step(const unsigned int t_mdl) const { return a.on_step(t_mdl) || b.on_step(t_mdl); }

  /**
   * @brief Run at the start of a step for combination of 2 observers.
   *
   * Each observer is run sequentially.
   *
   * @param t_mdl The unsigned int parameter.
   * @param d_gbxs The view of gridboxes in device memory.
   * @param totsupers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    a.at_start_step(t_mdl, d_gbxs, totsupers);
    b.at_start_step(t_mdl, d_gbxs, totsupers);
  }

  SDMMonitor auto get_monitor_of_sdm_processes() const {
    return NullSDMMonitor{};
  }  // TODO(CB) decide how to combine monitors
};

/**
 * @brief Overloaded operator >> to combine two Observers.
 *
 * @param obs1 First Observer.
 * @param obs2 Second Observer.
 * @return CombinedObserver<Obs1, Obs2> Combined Observer.
 */
auto operator>>(const Observer auto obs1, const Observer auto obs2) {
  return CombinedObserver(obs1, obs2);
}

/**
 * @brief Structure NullObserver does nothing at all.
 *
 * NullObserver defined for completion of Observer's Monoid Set.
 *
 */
struct NullObserver {
  /**
   * @brief No operations before timestepping.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {}

  /**
   * @brief No perations after timestepping.
   */
  void after_timestepping() const {}

  /**
   * @brief Next observation time is largest possible value.
   *
   * @param t_mdl Unsigned int for current timestep.
   * @return The next observation time (maximum unsigned int).
   */
  unsigned int next_obs(const unsigned int t_mdl) const { return LIMITVALUES::uintmax; }

  /**
   * @brief Check if on step always returns false.
   *
   * Null observer is never on_step.
   *
   * @param t_mdl The unsigned int parameter.
   * @return bool, always false.
   */
  bool on_step(const unsigned int t_mdl) const { return false; }

  /**
   * @brief No operations at the start of a step.
   *
   * @param t_mdl The unsigned int for the current timestep.
   * @param d_gbxs The view of gridboxes in device memory.
   * @param totsupers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {}

  SDMMonitor auto get_monitor_of_sdm_processes() const { return NullSDMMonitor{}; }
};

#endif  // LIBS_OBSERVERS_OBSERVERS_HPP_
