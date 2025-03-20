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
 * Last Modified: Thursday 6th March 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions for Kokkos compatibile thread-safe versions of C++ Fisher-Yates
 * serial shuffling algorithm, and of MergeShuffle parallelsised shuffling
 * algorithm. MergeShuffle comes from "A Very Fast, Parallel Random Permutation
 * Algorithm", Axel Bacher, Olivier Bodini, Alexandros Hollender, and Jérémie
 * Lumbroso, August 14, 2015. See also their code repository:
 * https://github.com/axel-bacher/mergeshuffle)
 */

#include "./shuffle.hpp"

/*
 * wrapper function to make it easier to chaneg shuffling algorithm (e.g. to
 * switch to Fisher-Yates for debugging)
 */
KOKKOS_FUNCTION viewd_supers shuffle_supers(const TeamMember& team_member,
                                            const viewd_supers supers,
                                            const GenRandomPool genpool) {
  if (team_member.team_size() > 1) {
    return merge_shuffle_supers(team_member, supers, genpool);  // parallelised shuffle
  } else {
    return fisher_yates_shuffle_supers(team_member, supers, genpool);  // default serial shuffle
  }
}

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
KOKKOS_FUNCTION viewd_supers fisher_yates_shuffle_supers(const TeamMember& team_member,
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

/*
 * C++ and Kokkos compatible version of ```shuffle(unsigned int *t, unsigned int
 * n)``` from "A Very Fast, Parallel Random Permutation Algorithm", Axel Bacher,
 * Olivier Bodini, Alexandros Hollender, and Jérémie Lumbroso, August 14, 2015.
 * see: https://github.com/axel-bacher/mergeshuffle/blob/master/merge_omp.c
 */
KOKKOS_FUNCTION viewd_supers merge_shuffle_supers(const TeamMember& team_member,
                                                  const viewd_supers supers,
                                                  const GenRandomPool genpool) {
  namespace KE = Kokkos::Experimental;

  constexpr unsigned int cutoff = 1024;  // below, c is largest integer < log_2(nn / cutoff)
  /**< determines smallest number of superdroplets for which merge sort
   * resorts to fisher-yates */

  const size_t nn = supers.extent(0);  // total number of superdroplets to shuffle
  unsigned int c = 0;                  // length of blocks for fisher-yates = (n/2^c) +/- 1
  while ((nn >> c) > cutoff) c++;      // c is largest number such that (n/2^c > cutoff)
  unsigned int q = 1 << c;             // number of blocks, q = 2^c

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, q), [=](const unsigned int i) {
    const size_t j = (nn * i) >> c;        // (nn+i)/2^c
    const size_t k = (nn * (i + 1)) >> c;  // (nn+nn*i)/2^c

    const auto first = KE::begin(supers);
    const auto block_first = first + j;  // distance from 1st to j'th element
    const auto block_dist = KE::distance(block_first, first + k - 1);  // distance from 1st to k'th

    URBG<ExecSpace> urbg{genpool.get_state()};
    fisher_yates_shuffle(urbg, block_first, block_dist);
    genpool.free_state(urbg.gen);
  });
  team_member.team_barrier();  // synchronise threads

  for (unsigned int p = 1; p < q; p += p) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, q / (2 * p)),
                         [=](const unsigned int _i) {
                           const auto i = _i * 2 * p;
                           const size_t j = (nn * i) >> c;
                           const size_t k = (nn * (i + p)) >> c;
                           const size_t l = (nn * (i + 2 * p)) >> c;
                           URBG<ExecSpace> urbg{genpool.get_state()};
                           merge_blocks(urbg, supers, j, k, l);
                           genpool.free_state(urbg.gen);
                         });
  }

  team_member.team_barrier();  // synchronise threads

  return supers;
}
