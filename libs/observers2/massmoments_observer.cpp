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
 * Last Modified: Friday 29th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality to calculate mass moments of droplet distribution in each gridbox in parallel
 */

#include "./massmoments_observer.hpp"

/* calculated 0th, 1st and 2nd moments of the (real)
droplet mass distribution in each gridbox, i.e. 0th, 3rd and 6th
moments of the droplet radius distribution for each gridbox.
Calculation is done for all gridboxes in parallel.
Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}  */
KOKKOS_FUNCTION
void calculate_massmoments(const TeamMember &team_member, const int ii,
                           const subviewd_constsupers supers,
                           Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                           Buffer<float>::mirrorviewd_buffer d_mom1,
                           Buffer<float>::mirrorviewd_buffer d_mom2) {
  const size_t nsupers(supers.extent(0));
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team_member, nsupers),
      KOKKOS_LAMBDA(const size_t kk, uint32_t &m0, float &m1, float &m2) {
        const auto xi = static_cast<double>(
            supers(kk).get_xi());  // cast multiplicity from unsigned int to double
        const auto mass = supers(kk).mass();
        m0 += static_cast<uint32_t>(supers(kk).get_xi());
        m1 += static_cast<float>(xi * mass);
        m2 += static_cast<float>(xi * mass * mass);
      },
      d_mom0(ii), d_mom1(ii), d_mom2(ii));  // {0th, 1st, 2nd} mass moments
}

/* calculated 0th, 1st and 2nd moment of the (real)
droplet mass distribution, i.e. 0th, 3rd and 6th
moment of the droplet radius distribution for one gridbox.
For all gridboxes in parallel */
void calculate_massmoments(const viewd_constgbx d_gbxs, Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                           Buffer<float>::mirrorviewd_buffer d_mom1,
                           Buffer<float>::mirrorviewd_buffer d_mom2) {
  const size_t ngbxs(d_gbxs.extent(0));
  Kokkos::parallel_for(
      "calculate_massmoments", TeamPolicy(ngbxs, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const TeamMember &team_member) {
        const int ii = team_member.league_rank();

        auto supers(d_gbxs(ii).supersingbx.readonly());
        calculate_massmoments(team_member, ii, supers, d_mom0, d_mom1, d_mom2);
      });
}

/* calculated 0th, 1st and 2nd moments of the (real)
raindroplet mass distribution for one gridbox, i.e. 0th, 3rd and 6th
moments of the raindroplet radius distribution for one gridbox.
A raindrop is droplet with a radius >= rlim = 40microns.
Kokkos::parallel_reduce([...]) is equivalent in serial to:
for (size_t kk(0); kk < supers.extent(0); ++kk){[...]} */
KOKKOS_FUNCTION
void calculate_massmoments_raindrops(const TeamMember &team_member, const int ii,
                                     const subviewd_constsupers supers,
                                     Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                                     Buffer<float>::mirrorviewd_buffer d_mom1,
                                     Buffer<float>::mirrorviewd_buffer d_mom2) {
  const size_t nsupers(supers.extent(0));
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team_member, nsupers),
      KOKKOS_LAMBDA(const size_t kk, uint32_t &m0, float &m1, float &m2) {
        const auto xi = static_cast<double>(
            supers(kk).get_xi());  // cast multiplicity from unsigned int to double
        const auto mass = supers(kk).mass();
        m0 += static_cast<uint32_t>(supers(kk).get_xi());
        m1 += static_cast<float>(xi * mass);
        m2 += static_cast<float>(xi * mass * mass);
      },
      d_mom0(ii), d_mom1(ii), d_mom2(ii));  // {0th, 1st, 2nd} mass moments
}

/* calculated 0th, 1st and 2nd moments of the (real)
raindroplet mass distribution in each gridbox, i.e. 0th, 3rd and 6th
moments of the raindroplet radius distribution for each gridbox.
A raindrop is droplet with a radius >= rlim = 40microns.
Calculation is done for all gridboxes in parallel.
Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}  */
void calculate_massmoments_raindrops(const viewd_constgbx d_gbxs,
                                     Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                                     Buffer<float>::mirrorviewd_buffer d_mom1,
                                     Buffer<float>::mirrorviewd_buffer d_mom2) {
  const size_t ngbxs(d_gbxs.extent(0));
  Kokkos::parallel_for(
      "calculate_massmoments", TeamPolicy(ngbxs, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const TeamMember &team_member) {
        const int ii = team_member.league_rank();

        auto supers(d_gbxs(ii).supersingbx.readonly());
        calculate_massmoments_raindrops(team_member, ii, supers, d_mom0, d_mom1, d_mom2);
      });
}
