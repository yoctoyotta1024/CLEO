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
#include <Kokkos_StdAlgorithms.hpp>

#include "./superdrop.hpp"
#include "./kokkosaliases_sd.hpp"

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
  Kokkos::Random_XorShift64<DeviceType> gen;

  KOKKOS_INLINE_FUNCTION
  uint64_t operator()(const uint64_t start,
                      const uint64_t end)
  /* draws a random 64 bit unsigned int from
  uniform distribution in the range [start, end] */
  {
    return gen.urand(start, end); // unsigned int rand
  }

  KOKKOS_INLINE_FUNCTION
  double drand(const double start,
               const double end)
  /* draws a random number (double) from uniform
  distribution in the range [0.0, 1.0] */
  {
    return gen.drand(start, end); // double rand
  }
};


KOKKOS_INLINE_FUNCTION
void device_swap(Superdrop& a, Superdrop& b)
/* swaps the values of the superdroplets a and b 
like C++98 std::swap except function works
on device as well as host. Note: Involves a copy
construction and two assignment operations
=> not efficient way of swapping the contents if 
Superdrop class stores large quantities of data */
{
  Superdrop c(a);
  a=b;
  b=c;
}

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
viewd_supers shuffle_supers(const viewd_supers supers,
                            URBG<DeviceType> urbg)
{
  namespace KE = Kokkos::Experimental;

  const auto first = KE::begin(supers);
  const auto dist = KE::distance(first, KE::end(supers) - 1); // distance to last element from first

  for (auto iter(dist); iter > 0; --iter)
  {
    const auto randiter = urbg(0, iter); // random uint64_t equidistributed between [0, i]
    device_swap(*(first + iter), *(first + randiter));
  }

  return supers;
}

KOKKOS_INLINE_FUNCTION viewd_supers
one_shuffle_supers(const TeamMember &team_member,
                   const viewd_supers supers,
                   GenRandomPool genpool)
/* Use only 1 member of team to randomly shuffle order
of superdroplet objects in 'supers'. Then synchronise
team and return view of shuffled supers. Uses 
to get team / thread safe random number generator */
{
  Kokkos::single(Kokkos::PerTeam(team_member), [&]()
                 {
                    URBG<ExecSpace> urbg{genpool.get_state()};
                    shuffle_supers(supers, urbg);
                    genpool.free_state(urbg.gen); });

  team_member.team_barrier(); // synchronize threads

  return supers;
}

#endif // URBG_HPP