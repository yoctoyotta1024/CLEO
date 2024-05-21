/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_condensation.hpp
 * Project: sdmmonitor
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
 * struct to monitor condensation SDM microphysical process
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_MONITOR_CONDENSATION_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_MONITOR_CONDENSATION_HPP_

#include <Kokkos_Core.hpp>

struct MonitorCondensation {
  Kokkos::View<double[1]> condrate;  // TODO(CB) monitor condensation properly

  MonitorCondensation() : condrate("condrate") { Kokkos::deep_copy(condrate, 0.0); }

  void monitor_microphysics() const;
};

#endif  //  LIBS_OBSERVERS_SDMMONITOR_MONITOR_CONDENSATION_HPP_
