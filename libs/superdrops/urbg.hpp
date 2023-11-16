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
satisfy requirements of C++11 UniformRandomBitGenerator
bject for a 32 bit unsigned int. Useful e.g. so that
gen's urand() function can be used in std::shuffle
to generate random pairs of superdroplets
during collision process */
{
  using result_type = uint64_t;
  Kokkos::Random_XorShift64<DeviceType> gen;
  
  static constexpr result_type min()
  {
    return LIMITVALUES::uint64tmin;
  }
  static constexpr result_type max()
  /* is equivalent to return
  Kokkos::Random_XorShift64<DeviceType>::MAX_URAND; */
  {
    return LIMITVALUES::uint64tmax;
  }

  result_type operator()()
  /* draws a random number from uniform
  distribution in the range [0,MAX_URAND] */
  {
    return gen.urand();
  }

  result_type operator()(const uint64_t range)
  /* draws a random number from uniform
  distribution in the range [0, range] */
  {
    return gen.urand(range);
  }

  result_type operator()(const uint64_t start,
                         const uint64_t end)
  /* draws a random number from uniform
  distribution in the range [start, end] */
  {
    return gen.urand(start, end);
  }
};

#endif // URBG_HPP