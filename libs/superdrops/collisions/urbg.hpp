/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: urbg.hpp
 * Project: collisions
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Struct (for Kokkos compatibility) to generate random numbers for SDM (e.g. to randomly shuffle
 * super-droplet's vector) based on C++11 standard UniformRandomBitGenerator (URBG). File also
 * contains Kokkos compatibile thread-safe versions of C++ shuffling algorithms.
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_URBG_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_URBG_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <cstdint>

/**
 * @brief Struct wrapping Kokkos random number generator.
 *
 * Generates random numbers in the range [start, end). Result equivalent to
 * std::uniform_int_distribution with parameters [a, b) = [start, end), where g = C++11
 * UniformRandomBitGenerator (URBG). Useful e.g. for using urand(start, end) function of Kokkos
 * random number generator 'gen' to generate random numbers for shuffling super-droplets array by
 * swapping elements in range [start, end), e.g, for linear sampling of super-droplet pairs in
 * SDM collision algorithm.
 *
 * @tparam DeviceType The Kokkos device type.
 */
template <class DeviceType>
struct URBG {
  Kokkos::Random_XorShift64<DeviceType> gen; /**< Kokkos random number generator */

  /**
   * @brief Draws a random 64-bit unsigned integer (uint64_t) from a uniform distribution in the
   * range [start, end).
   *
   * includes start, excludes end value.
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
   * @brief Draws a random number (double) from a uniform distribution in the range [start, end).
   *
   * includes start, excludes end value.
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

#endif  // LIBS_SUPERDROPS_COLLISIONS_URBG_HPP_
