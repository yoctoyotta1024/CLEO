/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: sdmmonitor.hpp
 * Project: observers
 * Created Date: Wednesday 8th May 2024
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
 * concept and structs used by observers to monitor various SDM processes
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_HPP_

/**
 * @brief Concept of SDMmonitor to monitor various SDM processes.
 *
 * @tparam SDMMo Type that satisfies the SDMMonitor concept.
 */
template <typename SDMMo>
concept SDMMonitor = requires(SDMMo mo) {
  { mo.monitor_microphysics() } -> std::same_as<void>;
};

struct NullSDMMonitor {
  double WIP; /**< work in progress TODO(CB) Note: must be GPU compatible */

  NullSDMMonitor() = default;   // Kokkos requirement for a (dual)View
  ~NullSDMMonitor() = default;  // Kokkos requirement for a (dual)View

  void monitor_microphysics() const {}
};

#endif  // LIBS_OBSERVERS_SDMMONITOR_HPP_
