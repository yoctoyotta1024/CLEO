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
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Implementaiton file for functions for Kokkos compatibile thread-safe
 * versions of C++ Fisher-Yates serial shuffling algorithm.
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../kokkosaliases_sd.hpp"
#include "../superdrop.hpp"
#include "./urbg.hpp"

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
                                            const viewd_supers supers, const GenRandomPool genpool);
/**
 * @brief Swaps the values of two super-droplets.
 *
 * Equivalent to C++98 std::swap but works on device as well as host (gpu compatible).
 *
 * _Note:_ Involves a copy construction and two assignment operations, which may not be
 * efficient if Superdrop class stores large quantities of data.
 *
 * @param a The first super-droplet.
 * @param b The second super-droplet.
 */
KOKKOS_INLINE_FUNCTION void device_swap(Superdrop& a, Superdrop& b) {
  Superdrop c(a);
  a = b;
  b = c;
}

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
 * @param urbg The random number generator.
 * @param first iterator/pointer to first element in supers to shuffle
 * @param dist number of elements (including first) to shuffle
 * @return The shuffled view of super-droplets.
 */
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void fisher_yates_shuffle(URBG<DeviceType> urbg, const auto first,
                                                 const auto dist) {
  for (auto iter(dist); iter > 0; --iter) {
    const auto randiter = urbg(0, iter + 1);  // random uint64_t equidistributed between [0, i]
    device_swap(*(first + iter), *(first + randiter));
  }
}

#endif  // LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_
