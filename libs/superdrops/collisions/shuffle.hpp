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
 * Functions for Kokkos compatibile thread-safe versions of C++ Fisher-Yates serial shuffling
 * algorithm, and of MergeShuffle parallelsised shuffling algorithm. MergeShuffle comes from
 * "A Very Fast, Parallel Random Permutation Algorithm", Axel Bacher, Olivier Bodini,
 * Alexandros Hollender, and Jérémie Lumbroso, August 14, 2015. See also their code
 * repository: https://github.com/axel-bacher/mergeshuffle)
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "../kokkosaliases_sd.hpp"
#include "../superdrop.hpp"
#include "./urbg.hpp"

struct FisherYatesShuffle {
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
  KOKKOS_INLINE_FUNCTION void device_swap(Superdrop& a, Superdrop& b) const {
    Superdrop c(a);
    a = b;
    b = c;
  }

  /**
   * @brief Shuffles the order of super-droplets in a view.
   *
   * Randomly shuffles the order of super-droplets using the
   * URBG (Uniform Random Bit Generator) struct (on device).
   *
   * @tparam DeviceType The Kokkos device type.
   * @param supers The view of super-droplets to shuffle.
   * @param urbg The random number generator.
   * @return The shuffled view of super-droplets.
   */
  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION viewd_supers shuffle_supers(const viewd_supers supers,
                                                     URBG<DeviceType> urbg) const {
    namespace KE = Kokkos::Experimental;

    const auto first = KE::begin(supers);
    const auto dist =
        KE::distance(first, KE::end(supers) - 1);  // distance to last element from 1st

    for (auto iter(dist); iter > 0; --iter) {
      const auto randiter = urbg(0, iter);  // random uint64_t equidistributed between [0, i]
      device_swap(*(first + iter), *(first + randiter));
    }

    return supers;
  }

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
  KOKKOS_INLINE_FUNCTION viewd_supers shuffle_supers(const TeamMember& team_member,
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
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_SHUFFLE_HPP_
