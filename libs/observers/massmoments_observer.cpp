/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: massmoments_observer.cpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality to calculate mass moments of (rain)droplet
 * distribution in each gridbox in parallel
 */

#include "./massmoments_observer.hpp"

/**
 * @brief Performs calculation of 0th, 1st, and 2nd moments of the (real)
 * droplet mass distribution for a single gridbox through reduction over super-droplets.
 *
 * This operator is a functor to perform the calculation of the 0th, 1st, and 2nd moments
 * of the droplet mass distribution in a gridbox (i.e. 0th, 3rd, and 6th moments of the
 * droplet radius distribution) within a Kokkos::parallel_reduce range policy
 * loop over superdroplets within a team policy loop over gridboxes.
 *
 * Kokkos::parallel_reduce([...]) is equivalent in serial to sum over result of:
 * for (size_t kk(0); kk < supers.extent(0); ++kk){[...]}.
 *
 * _Note:_ conversion from 8 to 4-byte precision for all mass moments: mom0 from size_t
 * (architecture dependent usually long unsigned int = 8 bytes) to 8 byte unsigned integer, and
 * mom1 and mom2 from double (8 bytes) to float (4 bytes).
 *
 * @param team_member The Kokkos team member.
 * @param supers The view of super-droplets for a gridbox (on device).
 * @param mom0 Reference to where to place value of 0th mass moment.
 * @param mom1 Reference to where to place value of 1st mass moment.
 * @param mom2 Reference to where to place value of 2nd mass moment.
 */
KOKKOS_FUNCTION
void calculate_massmoments(const TeamMember &team_member, const viewd_constsupers supers,
                           uint64_t &mom0, float &mom1, float &mom2) {
  const size_t nsupers(supers.extent(0));

  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team_member, nsupers),
      KOKKOS_LAMBDA(const size_t kk, uint64_t &m0, float &m1, float &m2) {
        const auto &drop(supers(kk));

        assert((drop.get_xi() < LIMITVALUES::uint64_t_max) &&
               "superdroplet mulitiplicy too large to represent with 4 byte unsigned integer");
        m0 += static_cast<uint64_t>(drop.get_xi());

        const auto mass = drop.mass();
        const auto xi = static_cast<double>(drop.get_xi());  // cast multiplicity to double
        m1 += static_cast<float>(xi * mass);
        m2 += static_cast<float>(xi * mass * mass);
      },
      mom0, mom1, mom2);  // {0th, 1st, 2nd} mass moments
}

/**
 * @brief Performs calculation of 0th, 1st, and 2nd moments of the (real) raindroplet mass
 * distribution in each gridbox.
 *
 * This operator is a functor to perform the calculation of the 0th, 1st, and 2nd moments
 * of the raindroplet mass distribution in each gridbox (i.e. 0th, 3rd, and 6th moments of the
 * droplet radius distribution) within a Kokkos::parallel_for range policy
 * loop over superdroplets within a team policy loop over gridboxes.
 *
 * A raindroplet is a droplet with a radius >= rlim = 40microns.
 *
 * Kokkos::parallel_reduce([...]) is equivalent in serial to sum over result of:
 * for (size_t kk(0); kk < supers.extent(0); ++kk){[...]}.
 *
 * _Note:_ conversion from 8 to 4-byte precision for all mass moments: mom0 from size_t
 * (architecture dependent usually long unsigned int = 8 bytes) to 8 byte unsigned integer, and
 * mom1 and mom2 from double (8 bytes) to float (4 bytes).
 *
 * @param team_member The Kokkos team member.
 * @param supers The view of super-droplets for a gridbox (on device).
 * @param mom0 Reference to where to place value of 0th mass moment.
 * @param mom1 Reference to where to place value of 1st mass moment.
 * @param mom2 Reference to where to place value of 2nd mass moment.
 */
KOKKOS_FUNCTION
void calculate_rainmassmoments(const TeamMember &team_member, const viewd_constsupers supers,
                               uint64_t &mom0, float &mom1, float &mom2) {
  constexpr double rlim(40e-6 / dlc::R0);  // dimless minimum radius of raindrop

  const size_t nsupers(supers.extent(0));
  Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team_member, nsupers),
      KOKKOS_LAMBDA(const size_t kk, uint64_t &m0, float &m1, float &m2) {
        const auto &drop(supers(kk));
        const auto binary = bool{drop.get_radius() >= rlim};  // 1 if droplet is raindrop, else 0

        assert((drop.get_xi() < LIMITVALUES::uint64_t_max) &&
               "superdroplet mulitiplicy too large to represent with 4 byte unsigned integer");
        m0 += static_cast<uint64_t>(drop.get_xi() * binary);

        const auto mass = drop.mass();
        const auto xi = static_cast<double>(drop.get_xi());  // cast multiplicity to double
        m1 += static_cast<float>(xi * mass * binary);
        m2 += static_cast<float>(xi * mass * mass * binary);
      },
      mom0, mom1, mom2);  // {0th, 1st, 2nd} mass moments
}
