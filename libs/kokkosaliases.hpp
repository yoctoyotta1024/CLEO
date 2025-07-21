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
#include <Kokkos_ScatterView.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include <memory>

#include "gridboxes/gridbox.hpp"
#include "superdrops/kokkosaliases_sd.hpp"
#include "superdrops/superdrop.hpp"

/* Gridboxes */
using dualview_gbx = Kokkos::DualView<Gridbox *>;            /**< Dual view of gridboxes. */
using dualview_constgbx = Kokkos::DualView<const Gridbox *>; /**< Dual view of const gridboxe. */

using viewh_gbx = dualview_gbx::t_host;           /**< View in host memory of gridboxes. */
using viewh_constgbx = dualview_constgbx::t_host; /**< view in host memory of const gridboxes. */

using viewd_gbx = dualview_gbx::t_dev;           /**< View in device memory of gridboxes. */
using viewd_constgbx = dualview_constgbx::t_dev; /**< View in device memory of const gridboxes. */

/* Gridbox Maps */
using kokkos_pairmap = Kokkos::UnorderedMap<unsigned int, Kokkos::pair<double, double>, ExecSpace>;
/**< E.g. for map from unsigned int gbxindex to gridbox boundaries */
using kokkos_uintmap = Kokkos::UnorderedMap<unsigned int, unsigned int, ExecSpace>;
/**< E.g. for map from one unsigned int gbxindex to another */
using kokkos_dblmaph = Kokkos::UnorderedMap<unsigned int, double, HostSpace>;
/**< E.g. for map from unsigned int gbxindex to gridbox area/volume on host*/
using viewd_ndims = Kokkos::View<size_t[3]>;
/**< View in device memory for number of gridboxes in CartesianMaps. */

/* Sorting Superdrops */
using viewd_counts = Kokkos::View<size_t *>; /**< View in device memory for sorting superdroplets */
/**< Scatter view for abstracted use of atomics/duplicates when computing sums for viewd_counts */
using scatterviewd_counts = Kokkos::Experimental::ScatterView<size_t *>;

namespace KokkosCleoSettings {
constexpr auto team_size = Kokkos::AUTO();
/**< configurable number threads per team for hierarchical parallelism over superdroplets */
}  // namespace KokkosCleoSettings

#endif  // LIBS_KOKKOSALIASES_HPP_
