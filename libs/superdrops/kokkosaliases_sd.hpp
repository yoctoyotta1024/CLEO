/*
 * ----- CLEO -----
 * File: kokkosaliases_sd.hpp
 * Project: superdrops
 * Created Date: Saturday 14th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * aliases for Kokkos superdrop views
 */

#ifndef KOKKOSALIASES_SD_HPP
#define KOKKOSALIASES_SD_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Random.hpp>

#include "./superdrop.hpp"

/* Default Execution Space for Parallelism */
// using ExecSpace = Kokkos::DefaultExecutionSpace;
using ExecSpace = Kokkos::DefaultHostExecutionSpace;
using HostSpace = Kokkos::DefaultHostExecutionSpace;

/* Superdrop views and subviews */
using viewd_supers = Kokkos::View<Superdrop *, ExecSpace::memory_space>;            // view in execution space's memory of superdroplets 
using viewd_constsupers = Kokkos::View<const Superdrop *, ExecSpace::memory_space>; // view in execution space's memory of const superdroplets

using subviewd_supers = Kokkos::Subview<viewd_supers, Kokkos::pair<size_t, size_t>>; // subiew of supers (for instance in a gridbox)
using subviewd_constsupers = Kokkos::Subview<viewd_constsupers, Kokkos::pair<size_t, size_t>>; // const supers subview (for instance in a gridbox) 

using mirrorh_constsupers = subviewd_constsupers::HostMirror; // mirror view (copy) of subview of superdroplets on host memory

/* Random Number Generation */
using GenRandomPool = Kokkos::Random_XorShift64_Pool<ExecSpace>; // type for pool of thread safe random number generators

/* Nested Parallelism */
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
using TeamMember = TeamPolicy::member_type;

#endif // KOKKOSALIASES_SD_HPP