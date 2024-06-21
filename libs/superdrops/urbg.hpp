/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: urbg.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Struct (for Kokkos compatibility) to generate random numbers for SDM (e.g. to randomly shuffle
 * super-droplet's vector) based on C++11 standard UniformRandomBitGenerator (URBG). File also
 * contains Kokkos compatibile thread-safe versions of C++ shuffling algorithms.
 */

#ifndef LIBS_SUPERDROPS_URBG_HPP_
#define LIBS_SUPERDROPS_URBG_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "kokkosaliases_sd.hpp"
#include "superdrop.hpp"

/**
 * @brief Struct wrapping Kokkos random number generator.
 *
 * Generates random numbers in the range [start, end]. Result equivalent to
 * std::uniform_int_distribution with parameters [a, b] = [start, end], where g = C++11
 * UniformRandomBitGenerator (URBG). Useful e.g. for using urand(start, end) function of Kokkos
 * random number generator 'gen' to generate random numbers for shuffling super-droplets array by
 * swapping elements in range [start, end] (e.g, for linear sampling of super-droplet pairs in
 * SDM collision algorithm).
 *
 * @tparam DeviceType The Kokkos device type.
 */
template <class DeviceType>
struct URBG {
  Kokkos::Random_XorShift64<DeviceType> gen; /**< Kokkos random number generator */

  /**
   * @brief Draws a random 64-bit unsigned integer (uint64_t) from a uniform distribution in the
   * range [start, end].
   *
   * @param start The lower bound of the range.
   * @param end The upper bound of the range.
   * @return The random 8-byte unsigned integer.
   */
  KOKKOS_INLINE_FUNCTION
  uint64_t operator()(const uint64_t start, const uint64_t end) {
    return gen.urand(start, end);  // unsigned int rand
  }

  /**
   * @brief Draws a random number (double) from a uniform distribution in the range [start, end].
   *
   * @param start The lower bound of the range.
   * @param end The upper bound of the range.
   * @return The random double.
   */
  KOKKOS_INLINE_FUNCTION
  double drand(const double start, const double end) {
    return gen.drand(start, end);  // double rand
  }
};

/**
 * @brief Swaps the values of two super-droplets.
 *
 * Equivalent to C++98 std::swap but works on device as well as host (gpu compatible).
 *
 * _Note:_ Involves a copy construction and two assignment operations, which may not be efficient
 * if Superdrop class stores large quantities of data.
 *
 * @param a The first super-droplet.
 * @param b The second super-droplet.
 */
KOKKOS_INLINE_FUNCTION
void device_swap(Superdrop& a, Superdrop& b) {
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
                                                   URBG<DeviceType> urbg) {
  namespace KE = Kokkos::Experimental;

  const auto first = KE::begin(supers);
  const auto dist = KE::distance(first, KE::end(supers) - 1);  // distance to last element from 1st

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
 * in the 'supers' view. Then synchronizes the team and returns the view of shuffled super-droplets.
 * Uses Kokkos thread safe random number generator from pool.
 *
 * @param team_member The Kokkos team member.
 * @param supers The view of superdroplets to shuffle.
 * @param genpool The random number generator pool.
 * @return The shuffled view of superdroplets.
 */
KOKKOS_INLINE_FUNCTION viewd_supers one_shuffle_supers(const TeamMember& team_member,
                                                       const viewd_supers supers,
                                                       const GenRandomPool genpool) {
  Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
    URBG<ExecSpace> urbg{genpool.get_state()};
    shuffle_supers(supers, urbg);
    genpool.free_state(urbg.gen);
  });

  team_member.team_barrier();  // synchronise threads

  return supers;
}

#endif  // LIBS_SUPERDROPS_URBG_HPP_
