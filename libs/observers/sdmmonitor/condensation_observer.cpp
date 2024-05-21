/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: condensation_observer.cpp
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
 * functionality to monitor condensation SDM microphysical process
 */

#include "./condensation_observer.hpp"

/**
 * @brief Monitor condensation rate
 *
 * _Note:_ conversion of condensation rate from double precision (8 bytes double) to single
 * precision (4 bytes float) in output.
 *
 */
KOKKOS_FUNCTION
void MonitorCondensation::monitor_microphysics() const {
  const auto rate_dbl = double{5.0};
  const auto rate = static_cast<float>(rate_dbl);
  Kokkos::deep_copy(d_data, rate);
}
