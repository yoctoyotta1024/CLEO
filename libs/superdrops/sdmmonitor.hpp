/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: sdmmonitor.hpp
 * Project: superdrops
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 27th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept and structs used by observers to monitor various SDM processes
 */

#ifndef LIBS_SUPERDROPS_SDMMONITOR_HPP_
#define LIBS_SUPERDROPS_SDMMONITOR_HPP_

#include <concepts>

#include "./kokkosaliases_sd.hpp"

/**
 * @brief Concept of SDMmonitor to monitor various SDM processes.
 *
 * @tparam SDMMo Type that satisfies the SDMMonitor concept.
 */
template <typename SDMMo>
concept SDMMonitor =
    requires(SDMMo mo, const TeamMember &tm, const double d, const viewd_constsupers supers) {
      { mo.reset_monitor() } -> std::same_as<void>;
      { mo.monitor_microphysics(tm, d) } -> std::same_as<void>;
      { mo.monitor_microphysics(tm, supers) } -> std::same_as<void>;
      { mo.monitor_motion(tm, supers) } -> std::same_as<void>;
    };

/**
 * @brief Structure CombinedSDMMonitor represents a new monitor formed from combination of two
 * SDMMonitors 'a' and 'b'.
 *
 * @tparam SDMMo1 Type satisfying the SDMMonitor concept.
 * @tparam SDMMo2 Type satisfying the SDMMonitor concept.
 */
template <SDMMonitor SDMMo1, SDMMonitor SDMMo2>
struct CombinedSDMMonitor {
 private:
  SDMMo1 a; /**< First Monitor. */
  SDMMo2 b; /**< Second Monitor. */

 public:
  /**
   * @brief Construct a new CombinedSDMMonitor object.
   *
   * @param mo1 First Monitor.
   * @param mo2 Second Monitor.
   */
  CombinedSDMMonitor(const SDMMo1 mo1, const SDMMo2 mo2) : a(mo1), b(mo2) {}

  /**
   * @brief reset monitor for combination of 2 sdm monitors.
   *
   * Each monitor is reset sequentially.
   */
  void reset_monitor() const {
    a.reset_monitor();
    b.reset_monitor();
  }

  /**
   * @brief monitor microphysics for combination of 2 sdm monitors.
   *
   * Each monitor is run sequentially.
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember &tm, const double d) const {
    a.monitor_microphysics(tm, d);
    b.monitor_microphysics(tm, d);
  }

  /**
   * @brief monitor microphysics for combination of 2 sdm monitors.
   *
   * Each monitor is run sequentially.
   */
  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember &tm, const viewd_constsupers supers) const {
    a.monitor_microphysics(tm, supers);
    b.monitor_microphysics(tm, supers);
  }

  /**
   * @brief monitor motion for combination of 2 sdm monitors.
   *
   * Each monitor is run sequentially.
   */
  KOKKOS_FUNCTION
  void monitor_motion(const TeamMember &tm, const viewd_constsupers supers) const {
    a.monitor_microphysics(tm, supers);
    b.monitor_microphysics(tm, supers);
  }
};

/**
 * @brief Null monitor for SDM processes from observer.
 *
 * NullSDMMonitor does nothing
 */
struct NullSDMMonitor {
  void reset_monitor() const {}

  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember &team_member, const double d) const {}

  KOKKOS_FUNCTION
  void monitor_microphysics(const TeamMember &team_member, const viewd_constsupers supers) const {}

  KOKKOS_FUNCTION
  void monitor_motion(const TeamMember &team_member, const viewd_constsupers supers) const {}
};

#endif  //  LIBS_SUPERDROPS_SDMMONITOR_HPP_
