/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: runstats_observer.hpp
 * Project: observers
 * Created Date: Monday 8th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 25th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for making and ouputting
 * statistics related to runtime performance
 * e.g. of timestepping CLEO
 */

#ifndef LIBS_OBSERVERS_RUNSTATS_OBSERVER_HPP_
#define LIBS_OBSERVERS_RUNSTATS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/sdmmonitor.hpp"

/**
 * @struct RunStats
 * @brief Structure to hold runtime statistics.
 */
struct RunStats {
 private:
  Kokkos::Timer kokkostimer;

 public:
  double t0;      /**< Time of observer creation. */
  double t_start; /**< Time before timestepping run. */
  double t_end;   /**< Time at the end of timestepping. */

  /**
   * @brief Constructor for RunStats.
   */
  RunStats() : kokkostimer(Kokkos::Timer()), t0(0.0), t_start(0.0), t_end(0.0) {
    t0 = kokkostimer.seconds();
  }

  /**
   * @brief Returns time elapsed since kokkostimer was at time=t0.
   * @return Time elapsed since t0 in seconds.
   */
  double time_elapsed() const { return kokkostimer.seconds() - t0; }
};

/**
 * @class RunStatsObserver
 * @brief Class that satisfies the observer concept and makes and
 * outputs runtime performance statistics.
 */
class RunStatsObserver {
 private:
  unsigned int interval;                /**< Timestep between runtime observations. */
  std::shared_ptr<RunStats> stats;      /**< Pointer to runtime statistics. */
  std::filesystem::path stats_filename; /**< Filename to output runtime statistics. */

  /**
   * @brief Function to be executed at the start of each timestep.
   *
   * Plug function does nothing but exists in case of need to add functionality
   * at the start of a timestep.
   */
  void at_start_step() const {}

  /**
   * @brief Prints a summary of runtime statistics to the terminal window.
   */
  void print_summary() const;

  /**
   * @brief Writes out some of the runtime statistics to a file.
   *
   * Writes timing statistics out to a text file called stats_filename.
   */
  void write_to_file() const;

 public:
  /**
   * @brief Constructor for RunStatsObserver.
   * @param obsstep Model timestep between runstats observations.
   * @param stats_filename Filename to output run statistics.
   */
  RunStatsObserver(const unsigned int obsstep, const std::filesystem::path stats_filename)
      : interval(obsstep), stats(std::make_shared<RunStats>()), stats_filename(stats_filename) {}

  /**
   * @brief Records statistics before timestepping.
   *
   * Records time of t_start = time of call to before_timestepping.
   *
   * @param d_gbxs View of grid boxes.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    stats->t_start = stats->time_elapsed();
  }

  /**
   * @brief Records statistics after timestepping.
   *
   * Records time of t_end = time of call to after_timestepping.
   */
  void after_timestepping() const {
    stats->t_end = stats->time_elapsed();
    print_summary();
    write_to_file();
  }

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
   * @brief Executes statistics functionality at the start of each timestep.
   *
   * If timestep is on observation step, call the function to make a runstats observation.
   *
   * @param t_mdl Current model time.
   * @param d_gbxs View of grid boxes.
   * @param totsupers View of super grids.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    if (on_step(t_mdl)) {
      at_start_step();
    }
  }

  /**
   * @brief Get null monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer that does nothing
   */
  SDMMonitor auto get_sdmmonitor() const { return NullSDMMonitor{}; }
};

#endif  // LIBS_OBSERVERS_RUNSTATS_OBSERVER_HPP_
