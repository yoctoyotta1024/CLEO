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
 * Functions for Kokkos compatibile thread-safe versions of C++ Fisher-Yates
 * serial shuffling algorithm, and of MergeShuffle parallelsised shuffling
 * algorithm. MergeShuffle comes from "A Very Fast, Parallel Random Permutation
 * Algorithm", Axel Bacher, Olivier Bodini, Alexandros Hollender, and Jérémie
 * Lumbroso, August 14, 2015. See also their code repository:
 * https://github.com/axel-bacher/mergeshuffle)
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../kokkosaliases_sd.hpp"
#include "../superdrop.hpp"
#include "./urbg.hpp"

KOKKOS_FUNCTION viewd_supers merge_shuffle_supers(const TeamMember& team_member,
                                                  const viewd_supers supers,
                                                  const GenRandomPool genpool);

/*
 * wrapper function to make it easier to chaneg shuffling algorithm (e.g. to
 * switch to FisherYates for debugging)
 */
KOKKOS_INLINE_FUNCTION viewd_supers shuffle_supers(const TeamMember& team_member,
                                                   const viewd_supers supers,
                                                   const GenRandomPool genpool) {
  return merge_shuffle_supers(team_member, supers, genpool);
}

/**
 * @brief Swaps the values of two super-droplets.
 *
 * Equivalent to C++98 std::swap but works on device as well as host (gpu
 * compatible).
 *
 * _Note:_ Involves a copy construction and two assignment operations, which may
 * not be efficient if Superdrop class stores large quantities of data.
 *
 * @param a The first super-droplet.
 * @param b The second super-droplet.
 */
KOKKOS_INLINE_FUNCTION void device_swap(Superdrop& a, Superdrop& b) {
  Superdrop c(a);
  a = b;
  b = c;
}

/*
Thread-Safe Kokkos compatible version of Fisher-Yates shuffling algorithm for
super-droplets. Shuffling is done in serial (slow!) on single thread per team,
we maintain this code in case it's needed for easier debugging (not for
performance!)
*/
struct FisherYatesShuffle {
  /**
   * @brief Shuffles the order of super-droplets in a view.
   *
   * Randomly shuffles the order of super-droplets using the URBG
   * (Uniform Random Bit Generator) struct (on device). Supers included in
   * shuffle are from iterators in range [first, first+dist] (inclusive), e.g.
   * if first points to the 5th superdroplet and dist=2, then the 5th, 6th and
   * 7th superdroplets will be shuffled amongst each other.
   *
   * @tparam DeviceType The Kokkos device type.
   * @param supers The view of super-droplets to shuffle.
   * @param urbg The random number generator.
   * @param first iterator/pointer to first element in supers to shuffle
   * @param dist number of elements (including first) to shuffle
   * @return The shuffled view of super-droplets.
   */
  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION viewd_supers shuffle_supers(URBG<DeviceType> urbg,
                                                     const viewd_supers supers, const auto first,
                                                     const auto dist) const {
    for (auto iter(dist); iter > 0; --iter) {
      const auto randiter = urbg(0, iter + 1);  // random uint64_t equidistributed between [0, iter]
      device_swap(*(first + iter), *(first + randiter));
    }

    return supers;
  }

  /**
   * @brief Randomly shuffles the order of super-droplet objects in a view using
   * a single thread in a Kokkos team (i.e. a single team member).
   *
   * Uses only one member of a Kokkos team to randomly shuffle the order of
   * super-droplet objects in the 'supers' view using Kokkos compatible rewrite
   * of C++ standard library Fisher-Yates shuffling algorithm. Afterwards
   * synchronizes the team and returns the view of shuffled super-droplets. Uses
   * Kokkos thread safe random number generator from pool.
   *
   * @param team_member The Kokkos team member.
   * @param supers The view of superdroplets to shuffle.
   * @param genpool The random number generator pool.
   * @return The shuffled view of superdroplets.
   */
  KOKKOS_FUNCTION viewd_supers operator()(const TeamMember& team_member, const viewd_supers supers,
                                          const GenRandomPool genpool) const;
};

/*
 * C++ and Kokkos compatible version of ```merge(unsigned int *t, unsigned int
 * m, unsigned int n)``` from "A Very Fast, Parallel Random Permutation
 * Algorithm", Axel Bacher, Olivier Bodini, Alexandros Hollender, and Jérémie
 * Lumbroso, August 14, 2015. see:
 * https://github.com/axel-bacher/mergeshuffle/blob/master/merge.c
 */
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void merge_blocks(URBG<DeviceType> urbg, const viewd_supers supers,
                                         const size_t j, const size_t k, const size_t l) {
  namespace KE = Kokkos::Experimental;

  const auto first = KE::begin(supers);  // iterator to first superdrop (like C pointer in merge.c)
  auto u = KE::distance(first, first + j);  // initial start position of 1st block
  auto v = KE::distance(first, first + k);  // initial start position of 2nd block
  auto w = KE::distance(first, first + l);  // initial end position of 2nd block

  // take elements from two blocks until one block is exhausted
  while (true) {
    if (urbg.flip()) {
      if (v == w) {
        break;
      }
      device_swap(*(first + u), *(first + v));
      ++v;
    } else {
      if (u == v) {
        break;
      }
      // no swap nor increment here
    }
    ++u;
  }

  // finish merge with Fisher-Yates
  while (u < w) {
    const auto randiter = urbg(0, u + 1);  // random uint64_t equidistributed between [0, u]
    device_swap(*(first + randiter), *(first + u));
    ++u;
  }
}

#endif  // LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_
