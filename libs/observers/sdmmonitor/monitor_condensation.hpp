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

#include "zarr/buffer.hpp"

struct MonitorCondensation {
  Buffer<float>::mirrorviewd_buffer condrate;  // TODO(CB) monitor condensation properly

  MonitorCondensation() : condrate("condrate", 1) { Kokkos::deep_copy(condrate, 0.0); }

  /**
   * @brief Monitor condensation rate
   *
   * _Note:_ conversion of condensation rate from double precision (8 bytes double) to single
   * precision (4 bytes float) in output.
   *
   */
  KOKKOS_FUNCTION
  void monitor_microphysics() const;
};

#endif  //  LIBS_OBSERVERS_SDMMONITOR_MONITOR_CONDENSATION_HPP_
