/*
 * ----- CLEO -----
 * File: kokkosaliases.hpp
 * Project: libs
 * Created Date: Saturday 14th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 22nd October 2023
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

using dualview_gbx = Kokkos::DualView<Gridbox *>;            // dual view of gridboxes
using dualview_constgbx = Kokkos::DualView<const Gridbox *>; // dual view of const gridboxes

using viewh_gbx = dualview_gbx::t_host;           // view in host memory of gridboxes
using viewh_constgbx = dualview_constgbx::t_host; // view in host memory of const gridboxes

using viewd_gbx = dualview_gbx::t_dev;           // view in device memory of gridboxes
using viewd_constgbx = dualview_constgbx::t_dev; // view in device memory of const gridboxes

using viewd_supers = Kokkos::View<Superdrop *>;            // should match that in kokkosaliases.hpp
using viewd_constsupers = Kokkos::View<const Superdrop *>; // should match that in kokkosaliases.hpp

using subviewd_supers = Kokkos::Subview<viewd_supers, Kokkos::pair<size_t, size_t>>; // should match that in gridbox.hpp
using subviewd_constsupers = Kokkos::Subview<viewd_constsupers, Kokkos::pair<size_t, size_t>>; // should match that in gridbox.hpp

using mirrorh_constsupers = subviewd_constsupers::HostMirror; // mirror view (copy) of subview of superdroplets on host memory

#endif // KOKKOSALIASES_HPP