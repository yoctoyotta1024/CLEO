/*
 * ----- CLEO -----
 * File: kokkosaliases.hpp
 * Project: libs
 * Created Date: Saturday 14th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 19th October 2023
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

#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"

using dualview_gbx = Kokkos::DualView<Gridbox *>; // dual view of gridboxes

using viewh_gbx = Kokkos::View<Gridbox *>;        // view in host memory of gridboxes
using viewh_constgbx = Kokkos::View<const Gridbox *>;        // view in host memory of gridboxes

using viewd_gbx = Kokkos::View<Gridbox *>;        // view in device memory of gridboxes
using viewd_constgbx = Kokkos::View<const Gridbox *>;        // view in device memory of gridboxes

using viewd_supers = Kokkos::View<Superdrop *>;   // view in device memory of superdroplets
using viewd_constsupers = Kokkos::View<const Superdrop *>;   // view in device memory of superdroplets

using viewd_solute = Kokkos::View<               // view to solute properties
    std::shared_ptr<const SoluteProperties>[1]>; // (stored in device memory and acessed through shared pointer)

#endif // KOKKOSALIASES_HPP