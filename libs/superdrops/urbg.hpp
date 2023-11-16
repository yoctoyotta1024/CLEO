/*
 * ----- CLEO -----
 * File: urbg.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 16th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Struct (for Kokkos compatibility) to
 * generate random numbers for SDM (e.g. to
 * shuffle superdroplet's vector) based on
 * c++11 standard UniformRandomBitGenerator
 * (URBG)
 */

#ifndef URBG_HPP
#define URBG_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

template <class DeviceType>
struct URBG
/* struct wrapping Kokkos random number generator to
generate random 64 bit unsigned int in range [start, end].
Result is analogous to std::uniform_int_distribution with
params [a,b]=[start, end] and g = C++11 UniformRandomBitGenerator
is URBG operator called with (start, end) = (0, URAND_MAX).
Useful so that gen's urand(start, end) function can be used
to randomly shuffle a kokkos view by swapping elements 
in range [start, end] e.g. to generate random pairs of
superdroplets during collision process */
{
  using result_type = uint64_t;
  Kokkos::Random_XorShift64<DeviceType> gen;
  
  result_type operator()(const uint64_t start,
                         const uint64_t end)
  /* draws a random number from uniform
  distribution in the range [start, end] */
  {
    return gen.urand(start, end);
  }
};

template <class DeviceType>
viewd_supers shuffle_supers(const viewd_supers supers, URBG<DeviceType> urbg)
{
  namespace KE = Kokkos::Experimental;

  const auto first = KE::begin(supers);
  const auto last = KE::end(supers);
  const auto diff = KE::distance(first, last - 1);

  for (auto i(diff); i > 0; --i)
  {
    const auto randit = urbg(0, i); // random int equidistributed between [0, i]
    KE::iter_swap(first + i, first + randit);
  }

  return supers;
}

#endif // URBG_HPP