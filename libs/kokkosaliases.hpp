/*
 * ----- CLEO -----
 * File: kokkosaliases.hpp
 * Project: libs
 * Created Date: Saturday 14th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 2nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * aliases for Kokkos views / dualviews / vectors
 */

#ifndef KOKKOSALIASES_HPP
#define KOKKOSALIASES_HPP

#include <memory>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Random.hpp>

#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/kokkosaliases_sd.hpp"

using dualview_gbx = Kokkos::DualView<Gridbox *>;            // dual view of gridboxes
using dualview_constgbx = Kokkos::DualView<const Gridbox *>; // dual view of const gridboxes

using viewh_gbx = dualview_gbx::t_host;           // view in host memory of gridboxes
using viewh_constgbx = dualview_constgbx::t_host; // view in host memory of const gridboxes

using viewd_gbx = dualview_gbx::t_dev;           // view in device memory of gridboxes
using viewd_constgbx = dualview_constgbx::t_dev; // view in device memory of const gridboxes

using GenRandomPool = Kokkos::Random_XorShift64_Pool<ExecSpace>; // type for pool of thread safe random number generators

using viewd_ndims = Kokkos::View<size_t[3]>; // view in device memory for number of gridboxes in CartesianMaps

#endif // KOKKOSALIASES_HPP