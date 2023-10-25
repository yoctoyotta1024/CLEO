/*
 * ----- CLEO -----
 * File: urbg.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
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
  using result_type = uint32_t;
  Kokkos::Random_XorShift64<DeviceType> gen;

  URBG(Kokkos::Random_XorShift64<DeviceType> gen) : gen(gen){};
  
  static constexpr result_type min()
  {
    return LIMITVALUES::uint32tmin;
  }
  static constexpr result_type max()
  /* is equivalent to return
  Kokkos::Random_XorShift64<DeviceType>::MAX_URAND; */
  {
    return LIMITVALUES::uint32tmax;
  }

  result_type operator()()
  {
    return gen.urand();
  }
};

#endif // URBG_HPP