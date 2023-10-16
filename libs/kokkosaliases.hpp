/*
 * ----- CLEO -----
 * File: kokkosaliases.hpp
 * Project: libs
 * Created Date: Saturday 14th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 16th October 2023
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

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "sdmdomain/gridbox.hpp"
#include "superdrops/superdrop.hpp"

using dualview_gbx = Kokkos::DualView<Gridbox *>; // dual view of gridboxes

using viewh_gbx = Kokkos::View<Gridbox *>;        // view in host memory of gridboxes
using viewh_constgbx = Kokkos::View<const Gridbox *>;        // view in host memory of gridboxes

using viewd_gbx = Kokkos::View<Gridbox *>;        // view in device memory of gridboxes
using viewd_constgbx = Kokkos::View<const Gridbox *>;        // view in device memory of gridboxes

using viewd_supers = Kokkos::View<Superdrop *>;   // view in device memory of superdroplets

#endif // KOKKOSALIASES_HPP