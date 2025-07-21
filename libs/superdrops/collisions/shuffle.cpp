/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: shuffle.cpp
 * Project: collisions
 * Created Date: Thursday 6th March 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Implementaiton file for functions for Kokkos compatibile thread-safe
 * versions of C++ Fisher-Yates serial shuffling algorithm.
 */

#include "./shuffle.hpp"

/**
 * @brief Randomly shuffles the order of super-droplet objects in a view
 * using Fisher-Yates algorithm
 *
 * Thread-Safe Kokkos compatible version of Fisher-Yates shuffling algorithm for
 * super-droplets. Shuffling is done in serial (slow!) on single thread per team
 *
 * Uses only one member of a Kokkos team to randomly shuffle the order of super-droplet objects
 * in the 'supers' view using Kokkos compatible rewrite of C++ standard library Fisher-Yates
 * shuffling algorithm. Afterwards synchronizes the team and returns the view of
 * shuffled super-droplets. Uses Kokkos thread safe random number generator from pool.
 *
 * @param team_member The Kokkos team member.
 * @param supers The view of superdroplets to shuffle.
 * @param genpool The random number generator pool.
 * @return The shuffled view of superdroplets.
 */
KOKKOS_FUNCTION viewd_supers shuffle_supers(const TeamMember& team_member,
                                            const viewd_supers supers,
                                            const GenRandomPool genpool) {
  namespace KE = Kokkos::Experimental;

  Kokkos::single(Kokkos::PerTeam(team_member), [=]() {
    const auto first = KE::begin(supers);
    const auto dist = KE::distance(first, KE::end(supers) - 1);  // distance to last elemnt from 1st
    URBG<ExecSpace> urbg{genpool.get_state()};
    fisher_yates_shuffle(urbg, first, dist);
    genpool.free_state(urbg.gen);
  });
  team_member.team_barrier();  // synchronise threads

  return supers;
}
