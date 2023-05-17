// Author: Clara Bayley
// File: randomgen.hpp
/* Struct that controls how to generate
random numbers for SDM (e.g. to shuffle
superdroplet's vector) */

#ifndef RANDOMGEN_HPP
#define RANDOMGEN_HPP

#include <limits>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

template <class DeviceType>
struct URBG
/* struct staisfying requirement to meet requirements
of C++11 UniformRandomBitGenerator Object */
{
  using result_type = uint32_t;
  Kokkos::Random_XorShift64<DeviceType> &gen;

  URBG(Kokkos::Random_XorShift64<DeviceType> &gen) : gen(gen) {};

  static constexpr result_type min()
  {
    return std::numeric_limits<result_type>::min();
  }
  static constexpr result_type max()
  /* is equivalent to return 
  Kokkos::Random_XorShift64<DeviceType>::MAX_URAND; */
  {
    return std::numeric_limits<result_type>::max(); 
  }

  result_type operator()()
  {
    return gen.urand();
  }
};

#endif // RANDOMGEN_HPP