/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: shuffle.hpp
 * Project: collisions
 * Created Date: Thursday 6th March 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 6th March 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Implementaiton file for functions for Kokkos compatibile thread-safe versions of C++
 * Fisher-Yates serial shuffling algorithm, and of MergeShuffle parallelsised shuffling
 * algorithm. MergeShuffle comes from "A Very Fast, Parallel Random Permutation Algorithm",
 * Axel Bacher, Olivier Bodini, Alexandros Hollender, and Jérémie Lumbroso, August 14, 2015.
 * See also their code repository: https://github.com/axel-bacher/mergeshuffle)
 */

#include "./shuffle.hpp"

/**
 * @brief Randomly shuffles the order of super-droplet objects in a view using a single thread in
 * a Kokkos team (i.e. a single team member).
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
KOKKOS_FUNCTION viewd_supers FisherYatesShuffle::operator()(const TeamMember& team_member,
                                                            const viewd_supers supers,
                                                            const GenRandomPool genpool) const {
  Kokkos::single(Kokkos::PerTeam(team_member), [=, *this]() {
    URBG<ExecSpace> urbg{genpool.get_state()};
    shuffle_supers(supers, urbg);
    genpool.free_state(urbg.gen);
  });
  team_member.team_barrier();  // synchronise threads

  return supers;
}
