/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: streamout_observer.hpp
 * Project: observers
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 8th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Struct satisfies observer type and streams out live data to an output device
 * (e.g. computer screen) about gridboxes during every observation
 * at fixed 'interval' timesteps.
 */

#ifndef LIBS_OBSERVERS_STREAMOUT_OBSERVER_HPP_
#define LIBS_OBSERVERS_STREAMOUT_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./sdmmonitor.hpp"
#include "gridboxes/gridbox.hpp"

namespace dlc = dimless_constants;

/**
 * @struct StreamOutObserver
 * @brief Struct that satisfies the observer concept and streams out live data to an output device
 * (e.g. computer screen) about gridboxes during every observation at fixed 'interval' timesteps.
 */
struct StreamOutObserver {
 private:
  unsigned int interval; /**< Timestep between output events. */
  std::function<double(unsigned int)>
      step2realtime; /**< Function to convert model timesteps to real time. */

  /**
   * @brief Prints a statement about the state of grid boxes.
   *
   * This function prints out information about the state of gridboxes.
   * It extracts information from the 0th gridbox in the gridboxes' view and prints some information
   * e.g. temperature, pressure, specific humidity, and specific cloud water content.
   * Additionally, it prints the total number of superdroplets in the domain and the total number of
   * gridboxes.
   *
   * @param t_mdl Current model time.
   * @param d_gbxs View of the gridboxes on the device.
   */
  void streamout_statement(const unsigned int t_mdl, const viewd_constgbx d_gbxs) const;

 public:
  /**
   * @brief Constructor for StreamOutObserver.
   * @param obsstep Interval in model timesteps between observation events.
   * @param step2realtime Function to convert model timesteps to real time.
   */
  StreamOutObserver(const unsigned int obsstep,
                    const std::function<double(unsigned int)> step2realtime)
      : interval(obsstep), step2realtime(step2realtime) {}

  /**
   * @brief Function called before timestepping.
   * @param d_gbxs View of grid boxes.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes StreamOutObserver\n";
  }

  /**
   * @brief Function called after timestepping.
   */
  void after_timestepping() const {}

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
   * @brief Observe gridboxes at the start of each timestep.
   *
   * If timestep is on observation step, stream out statement about gridboxes to output device
   * e.g. a computer terminal.
   *
   * @param t_mdl Current model time.
   * @param d_gbxs View of grid boxes.
   * @param totsupers View of super grids.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    if (on_step(t_mdl)) {
      streamout_statement(t_mdl, d_gbxs);
    }
  }

  SDMMonitor get_monitor_of_sdm_processes() const { return SDMMonitor{}; }
};

#endif  // LIBS_OBSERVERS_STREAMOUT_OBSERVER_HPP_
