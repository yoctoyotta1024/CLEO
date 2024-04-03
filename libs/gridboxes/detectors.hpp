/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: detectors.hpp
 * Project: gridboxes
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 6th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for functions and structures related to detectors
 * which track data for output (e.g. of microphysical processes)
 * in gridboxes
 */


#ifndef LIBS_GRIDBOXES_DETECTORS_HPP_
#define LIBS_GRIDBOXES_DETECTORS_HPP_

#include <Kokkos_Core.hpp>

struct Detectors {
  Detectors() = default;   // Kokkos requirement for a (dual)View
  ~Detectors() = default;  // Kokkos requirement for a (dual)View
};

#endif  // LIBS_GRIDBOXES_DETECTORS_HPP_
