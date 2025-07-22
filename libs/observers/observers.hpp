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
#include "superdrops/sdmmonitor.hpp"

/**
 * @brief Concept Observer is all types that have functions for timestepping and
 * observations functions as constrained here.
 *
 * @tparam Obs The type that satisfies the Observer concept.
 */
template <typename Obs>
concept Observer = requires(Obs obs, unsigned int t, const viewd_constgbx d_gbxs,
                            const subviewd_constsupers d_supers) {
  { obs.next_obs(t) } -> std::convertible_to<unsigned int>;
  { obs.on_step(t) } -> std::same_as<bool>;
  { obs.before_timestepping(d_gbxs, d_supers) } -> std::same_as<void>;
  { obs.after_timestepping() } -> std::same_as<void>;
  { obs.at_start_step(t, d_gbxs, d_supers) } -> std::same_as<void>;
  { obs.get_sdmmonitor() };
};

/**
 * @brief Structure CombinedObserver represents a new observer formed from combination of two
 * Observers 'a' and 'b'.
 *
 * @tparam Obs1 Type satisfying the Observer concept.
 * @tparam Obs2 Type satisfying the Observer concept.
 */
template <Observer Obs1, Observer Obs2, SDMMonitor SDMMo>
struct CombinedObserver {
 private:
  Obs1 a;   /**< First Observer. */
  Obs2 b;   /**< Second Observer. */
  SDMMo mo; /**< Combination of First and Second Observers' SDMMonitors */

 public:
  /**
   * @brief Construct a new CombinedObserver object.
   *
   * @param obs1 First Observer.
   * @param obs2 Second Observer.
   * @param mo12 Combined Monitor from first and second observers.
   */
  CombinedObserver(const Obs1 obs1, const Obs2 obs2, const SDMMo mo12)
      : a(obs1), b(obs2), mo(mo12) {}

  /**
   * @brief Run before timestepping for combination of 2 observers.
   *
   * Each observer is run sequentially.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   */
  void before_timestepping(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    a.before_timestepping(d_gbxs, d_supers);
    b.before_timestepping(d_gbxs, d_supers);
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

    return Kokkos::min(t_a, t_b);
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
   * @param d_supers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const subviewd_constsupers d_supers) const {
    a.at_start_step(t_mdl, d_gbxs, d_supers);
    b.at_start_step(t_mdl, d_gbxs, d_supers);
  }

  /**
   * @brief Get monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer
   */
  SDMMonitor auto get_sdmmonitor() const { return mo; }
};

/**
 * @brief Return a new CombinedObserver object from combibation of 2 observers.
 *
 * @param obs1 First Observer.
 * @param obs2 Second Observer.
 */
Observer auto CombineObservers(const Observer auto obs1, const Observer auto obs2) {
  const SDMMonitor auto mo12 = CombinedSDMMonitor(obs1.get_sdmmonitor(), obs2.get_sdmmonitor());
  return CombinedObserver(obs1, obs2, mo12);
}

/**
 * @brief Overloaded operator >> to combine two Observers.
 *
 * @param obs1 First Observer.
 * @param obs2 Second Observer.
 * @return CombinedObserver<Obs1, Obs2> Combined Observer.
 */
auto operator>>(const Observer auto obs1, const Observer auto obs2) {
  return CombineObservers(obs1, obs2);
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
  void before_timestepping(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
  }

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
   * @param d_supers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const subviewd_constsupers d_supers) const {}

  /**
   * @brief Get null monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer that does nothing
   */
  SDMMonitor auto get_sdmmonitor() const { return NullSDMMonitor{}; }
};

#endif  // LIBS_OBSERVERS_OBSERVERS_HPP_
