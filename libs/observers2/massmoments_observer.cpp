/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: massmoments_observer.cpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 8th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality to calculate mass moments of (rain)droplet
 * distribution in each gridbox in parallel
 */

#include "./massmoments_observer.hpp"

/* Function performs calculation of 0th, 1st and 2nd moments of the (real) droplet mass distribution
in each gridbox, i.e. 0th, 3rd and 6th moments of the droplet radius distribution for each gridbox.
Calculation is done for all gridboxes in parallel.
Kokkos::parallel_reduce([...]) is equivalent in serial to:
for (size_t kk(0); kk < supers.extent(0); ++kk){[...]}.
Note conversion from 8 to 4byte precision for all mass moments: mom0 from size_t (architecture
dependent usually long unsigned int = 8 bytes) to single precision (uint32_t = 4 bytes), and mom1
and mom2 from double (8 bytes) to float (4 bytes) */
KOKKOS_FUNCTION
MassMomentsFunc::operator()(const TeamMember & team_member, const viewd_constgbx d_gbxs,
                            Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                            Buffer<float>::mirrorviewd_buffer d_mom1,
                            Buffer<float>::mirrorviewd_buffer d_mom2) {}

/* Function performs calculation of 0th, 1st and 2nd moments of the (real)
raindroplet mass distribution in each gridbox, i.e. 0th, 3rd and 6th moments of the raindroplet
radius distribution for each gridbox. A raindrop is droplet with a radius >= rlim = 40microns.
Calculation is done for all gridboxes in parallel.
Kokkos::parallel_reduce([...]) is equivalent in serial to:
for (size_t kk(0); kk < supers.extent(0); ++kk){[...]}.
Note conversion from 8 to 4byte precision for all mass moments: mom0 from size_t (architecture
dependent usually long unsigned int = 8 bytes) to single precision (uint32_t = 4 bytes), and mom1
and mom2 from double (8 bytes) to float (4 bytes) */
KOKKOS_FUNCTION
RaindropsMassMomentsFunc::operator()(const TeamMember & team_member, const viewd_constgbx d_gbxs,
                                     Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                                     Buffer<float>::mirrorviewd_buffer d_mom1,
                                     Buffer<float>::mirrorviewd_buffer d_mom2) {}
