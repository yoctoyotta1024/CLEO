/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: kokkosaliases.hpp
 * Project: libs
 * Created Date: Saturday 14th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 1st May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * aliases for Kokkos views / dualviews / vectors
 */

#ifndef LIBS_KOKKOSALIASES_HPP_
#define LIBS_KOKKOSALIASES_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include "gridboxes/gridbox.hpp"
#include "superdrops/kokkosaliases_sd.hpp"
#include "superdrops/superdrop.hpp"

/* Gridboxes */
using dualview_gbx = Kokkos::DualView<Gridbox *>;             // dualview of gridboxes
using dualview_constgbx = Kokkos::DualView<const Gridbox *>;  // dualview of const gridboxes

using viewh_gbx = dualview_gbx::t_host;            // view in host memory of gridboxes
using viewh_constgbx = dualview_constgbx::t_host;  // view in host memory of const gridboxes

using viewd_gbx = dualview_gbx::t_dev;            // view in device memory of gridboxes
using viewd_constgbx = dualview_constgbx::t_dev;  // view in device memory of const gridboxes

/* Gridbox Maps */
using viewd_ndims = Kokkos::View<size_t[3]>;
using kokkos_pairmap = Kokkos::UnorderedMap<unsigned int, Kokkos::pair<double, double>, ExecSpace>;
using kokkos_dblmap = Kokkos::UnorderedMap<unsigned int, double, ExecSpace>;
using kokkos_uintmap = Kokkos::UnorderedMap<unsigned int, unsigned int, ExecSpace>;

#endif  // LIBS_KOKKOSALIASES_HPP_
