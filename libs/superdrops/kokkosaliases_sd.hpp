/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: kokkosaliases_sd.hpp
 * Project: superdrops
 * Created Date: Saturday 14th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Aliases for Kokkos super-droplet views and parallelisation
 */

#ifndef LIBS_SUPERDROPS_KOKKOSALIASES_SD_HPP_
#define LIBS_SUPERDROPS_KOKKOSALIASES_SD_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_Random.hpp>

#include "superdrop.hpp"

/* Defines aliases for the (default Kokkos) execution spaces and memory spaces for parallelism. */
using ExecSpace =
    Kokkos::DefaultExecutionSpace; /**< (default) execution space for device parallelism */
using HostSpace =
    Kokkos::DefaultHostExecutionSpace; /**< (default) execution space for host parallelism */

/* Defines aliases for Kokkos views and subviews of super-droplets. */
using viewd_supers = Kokkos::View<Superdrop *>; /**< View in device memory of superdrops. */
using viewd_constsupers =
    Kokkos::View<const Superdrop *>; /**< View in device memory of const superdrops. */

using kkpair_size_t = Kokkos::pair<size_t, size_t>; /**< size_t pair (e.g. see supersingbx refs). */
using subviewd_supers =
    Kokkos::Subview<viewd_supers,
                    kkpair_size_t>; /**< Sub-View of supers (e.g. for instance in a gridbox). */
using subviewd_constsupers =
    Kokkos::Subview<viewd_constsupers,
                    kkpair_size_t>; /**< Const supers subview (e.g. for instance in a gridbox). */

/* Defines aliases for team policies and team members in heirarchal (i.e. nested) parallelism. */
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>; /**< Team policy in the execution space. */
using TeamMember = TeamPolicy::member_type;       /**< Member in device parallel execution team. */

using HostTeamPolicy = Kokkos::TeamPolicy<HostSpace>; /**< Team policy in the host space. */
using HostTeamMember = HostTeamPolicy::member_type;   /**< Member in host parallel execution team.*/

/**
 * @brief Thread-safe random number generation.
 *
 * Defines an alias for a pool of Kokkos random number generators in the execution space.
 */
using GenRandomPool = Kokkos::Random_XorShift64_Pool<ExecSpace>;

using viewd_constcoords = Kokkos::View<const double[3]>;
using viewd_coords = Kokkos::View<double[3]>;
/**< Helper type for view in device memory for superdroplet coords: (coord3, coord1, coord2). */

#endif  // LIBS_SUPERDROPS_KOKKOSALIASES_SD_HPP_
