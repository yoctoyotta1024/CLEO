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
 * Last Modified: Wednesday 22nd May 2024
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

/**
 * @brief Concept of SDMmonitor to monitor various SDM processes.
 *
 * @tparam SDMMo Type that satisfies the SDMMonitor concept.
 */
template <typename SDMMo>
concept SDMMonitor = requires(SDMMo mo) {
  { mo.monitor_microphysics() } -> std::same_as<void>;
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
   * @brief monitor microphysics for combination of 2 sdm monitors.
   *
   * Each monitor is run sequentially.
   */
  void monitor_microphysics() const {
    a.monitor_microphysics();
    b.monitor_microphysics();
  }
};

struct NullSDMMonitor {
  void monitor_microphysics() const {}
};

#endif  //  LIBS_SUPERDROPS_SDMMONITOR_HPP_
