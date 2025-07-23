/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_precipitation_observer.cpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to monitor precipitation during superdroplet motion
 */

#include "./monitor_precipitation_observer.hpp"

/**
 * @brief Parallel loop to fill d_data with zero value.
 */
void MonitorPrecipitation::reset_monitor() const {
  Kokkos::parallel_for(
      "reset_monitor", Kokkos::RangePolicy(0, d_data.extent(0)),
      KOKKOS_CLASS_LAMBDA(const size_t& jj) { d_data(jj) = 0.0; });
}
