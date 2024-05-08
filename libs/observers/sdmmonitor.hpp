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
 * Last Modified: Wednesday 8th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * TODO(CB) fill in struct for different monitors (and turn SDMMonitor into a concept?)
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_HPP_

struct SDMMonitor {
  double WIP; /**< work in progress TODO(CB) Note: must be GPU compatible */

  SDMMonitor() = default;   // Kokkos requirement for a (dual)View
  ~SDMMonitor() = default;  // Kokkos requirement for a (dual)View
};

#endif  // LIBS_OBSERVERS_SDMMONITOR_HPP_
