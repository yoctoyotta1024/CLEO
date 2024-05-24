/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: consttstep_observer.hpp
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
 * Concept and related structures for creating an observer which acts at the start of each step.
 */

#ifndef LIBS_OBSERVERS_CONSTTSTEP_OBSERVER_HPP_
#define LIBS_OBSERVERS_CONSTTSTEP_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "sdmmonitor/sdmmonitor.hpp"

/**
 * @brief Concept ObsFuncs for all types that can be called used by ConsttepObserver for
 * observation functions.
 *
 * Type in ConstTstepObserver obeying ObsFuncs makes it possible for ConstTstepObserver
 * to obey Observer concept.
 *
 * @tparam O Type that satisfies the ObsFuncs concept.
 */
template <typename OFs>
concept ObsFuncs = requires(OFs ofs, unsigned int t, const viewd_constgbx d_gbxs,
                            const viewd_constsupers totsupers) {
  { ofs.before_timestepping(d_gbxs) } -> std::same_as<void>;
  { ofs.after_timestepping() } -> std::same_as<void>;
  { ofs.at_start_step(t, d_gbxs, totsupers) } -> std::same_as<void>;
  { ofs.get_sdmmonitor() };
};

/**
 * @brief Structure ConstTstepObserver represents a type that satisfies the concept of an
 * observer with a constant timestep interval between observations at the start of each timestep.
 *
 * Struct can be used to create an observer with a constant timestep and with observation
 * functionality as determined by the 'do_obs' instance of the ObsFuncs type 'O'.
 *
 * @tparam O Type that satisfies the ObsFuncs concept.
 */
template <ObsFuncs O>
struct ConstTstepObserver {
 private:
  unsigned int interval; /**< interval between observations. */
  O do_obs;              /**< Observation functionality. */

 public:
  /**
   * @brief Construct a new ConstTstepObserver object.
   *
   * @param interval Timestep interval.
   * @param o Observer.
   */
  ConstTstepObserver(const unsigned int interval, const O o) : interval(interval), do_obs(o) {}

  /**
   * @brief Perform operations before timestepping.
   *
   * Calls `before_timestepping` function of `do_obs`.
   *
   * @param d_gbxs The view of gridboxes in device memory.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    do_obs.before_timestepping(d_gbxs);
  }

  /**
   * @brief Perform operations after timestepping.
   *
   * Calls `after_timestepping` function of `do_obs`.
   */
  void after_timestepping() const { do_obs.after_timestepping(); }

  /**
   * @brief Determine the next observation time.
   *
   * Calculates the next observation time based on the current model time and this observer's
   * constant timestep between observations, 'interval'.
   *
   * @param t_mdl The unsigned int parameter representing the current model timestep.
   * @return Unsigned int for the next observation timestep.
   */
  unsigned int next_obs(const unsigned int t_mdl) const {
    return ((t_mdl / interval) + 1) * interval;
  }

  /**
   * @brief Check if observer is "on step".
   *
   * Checks if the current model time is on an observation timestep.
   *
   * @param t_mdl The unsigned int parameter representing the current model timestep.
   * @return True if the current timestep is an observation timestep, false otherwise.
   */
  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  /**
   * @brief Perform operation at the start of a step if at appropriate interval.
   *
   * Calls `at_start_step` function of `do_obs` if the current model time is on
   * an observation timestep.
   *
   * @param t_mdl The unsigned int parameter representing the current model time.
   * @param d_gbxs The view of gridboxes in device memory.
   * @param totsupers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    if (on_step(t_mdl)) {
      do_obs.at_start_step(t_mdl, d_gbxs, totsupers);
    }
  }

  SDMMonitor auto get_sdmmonitor() const { return do_obs.get_sdmmonitor(); }
};

#endif  // LIBS_OBSERVERS_CONSTTSTEP_OBSERVER_HPP_
