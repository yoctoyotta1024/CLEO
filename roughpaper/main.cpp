/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 17th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * rough paper for checking small things
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <limits>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_NestedSort.hpp>


struct Superdrop
{
  size_t xi;
};

using viewd_supers = Kokkos::View<Superdrop *>;
using ExecSpace = Kokkos::DefaultExecutionSpace;

KOKKOS_INLINE_FUNCTION Kokkos::pair<Superdrop &, Superdrop &>
assign_drops(Superdrop &dropA, Superdrop &dropB)
/* compare dropA.xi with dropB.xi and return (non-const)
references to dropA and dropB in a pair {drop1, drop2}
such that drop1.xi is always >= drop2.xi */
{
  if (!(dropA.xi < dropB.xi))
  {
    return {dropA, dropB};
  }
  else
  {
    return {dropB, dropA};
  }
}

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

KOKKOS_INLINE_FUNCTION Kokkos::Subview<viewd_supers,
                                       Kokkos::pair<size_t, size_t>>
remove_null_supers(viewd_supers supers, const size_t nnull)
/* sort view of superdroplets by their sdgbxindexes
from lowest to highest sdgbxindex. Then set new subview
excluding supers with sdgbxindex = LIMITVALUES::uintmax */
{
  using TeamPol = Kokkos::TeamPolicy<ExecSpace>;

  struct SortComparator
  /* a precedes b if its sdgbxindex is smaller */
  {
    KOKKOS_INLINE_FUNCTION
    bool operator()(const Superdrop &a, const Superdrop &b) const
    {
      return (a.xi) < (b.xi);
    }
  };

  Kokkos::parallel_for(
      "sortingsupers_thread", TeamPol(1, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const TeamPol::member_type &t) {
        Kokkos::Experimental::sort_team(t, supers, SortComparator{});
      });

  const size_t nsupers(supers.extent(0) - nnull); // remove null supers from no. supers
  const Kokkos::pair<size_t, size_t> refs({0, nsupers});
  auto new_supers = Kokkos::subview(supers, refs);
  return new_supers;
}

int main(int argc, char *argv[])
{
  using GenRandomPool = Kokkos::Random_XorShift64_Pool<ExecSpace>; // type for pool of thread safe random number generators
 
  size_t nsupers(10);

  Kokkos::initialize(argc, argv);
  {
    GenRandomPool genpool(std::random_device{}());
    viewd_supers supers("supers", nsupers);

    auto h_supers = Kokkos::create_mirror_view(supers); // mirror of supers in case view is on device memory
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      h_supers(kk) = Superdrop{kk*10+kk};
    }
    Kokkos::deep_copy(supers, h_supers);

    {
      std::cout << " \n --- b4 ---\n ";
      for (size_t kk(0); kk < nsupers; ++kk)
      {
        std::cout << supers(kk).xi << ", ";
      }
      std::cout << " \n --- --- ---\n ";

      URBG<ExecSpace> urbg{genpool.get_state()}; // thread safe random number generator
      shuffle_supers(supers, urbg);

      std::cout << " \n --- l8r ---\n ";
      for (size_t kk(0); kk < nsupers; ++kk)
      {
        std::cout << supers(kk).xi << ", ";
      }
      std::cout << " \n --- --- ---\n ";

      supers = remove_null_supers(supers, 0);

      size_t new_nsupers(supers.extent(0));
      std::cout
          << " \n --- end ---\n ";
      for (size_t kk(0); kk < new_nsupers; ++kk)
      {
        std::cout << supers(kk).xi << ", ";
      }
      std::cout << " \n --- --- ---\n ";

      const double phi(urbg.drand(0.0, 1.0));
      std::cout << "\n phi: " << phi << "\n";
      genpool.free_state(urbg.gen);
    }
  }
  Kokkos::finalize();
}

// int main(int argc, char *argv[])
// {
//   std::vector<unsigned int> gbxidx = {0,10,20,30,40,50};
//   std::vector<double> bounds = {1.1,2.2,3.3,4.4,5.5,6.6};

//   const unsigned int idx(10);

//   auto it = std::find(gbxidx.begin(), gbxidx.end(), idx);
//   auto d = std::distance(gbxidx.begin(), it);

//   std::cout << "dis: " << d <<", "<< gbxidx.size() <<"\n";
//   if (d > (gbxidx.size()-1))
//   {
//     throw std::invalid_argument("idx not found in gbxidxs");
//   }
//   std::cout << gbxidx.at(d) << " -> " << bounds.at(d) << "\n";

//   return 0;
// }

// int main(int argc, char *argv[])
// {
//   Kokkos::initialize(argc, argv);
//   {

//     using stdp = Kokkos::pair<double, double>;

//     const unsigned int ngbxs = 3;

//     std::array<unsigned int, 3> keys = {2, 1, 0};
//     std::array<stdp, 3> vals = {stdp({-1.0, 1.0}),
//                                 stdp({-2.0, 2.0}),
//                                 stdp({-3.0, 3.0})};

//     Kokkos::UnorderedMap<unsigned int, stdp,
//                          Kokkos::DefaultExecutionSpace>
//         map4gbxs(ngbxs);

//     for (int i(0); i < ngbxs; ++i)
//     {
//       std::cout << "k: " << keys[i]
//                 << " -> ("
//                 << vals[i].first << " , "
//                 << vals[i].second << ")\n";

//       map4gbxs.insert(keys[i], vals[i]);
//     }

//     for (int i(0); i < ngbxs; ++i)
//     {
//       const unsigned int k(keys[i]);
//       const auto idx(map4gbxs.find(i));
//       std::cout << "idx: " << idx << "\n";
//       std::cout << "k: " << map4gbxs.key_at(idx)
//                 << " -> (" << map4gbxs.value_at(idx).first << ", "
//                 << map4gbxs.value_at(idx).second << ")\n";
//     }

//     std::cout << "\n---\n";
//     // assume umap is an existing Kokkos::UnorderedMap
//    Kokkos::parallel_for(
//         map4gbxs.capacity(), KOKKOS_LAMBDA(uint32_t i) {
//           if (map4gbxs.valid_at(i))
//           {
//             auto key = map4gbxs.key_at(i);
//             auto value = map4gbxs.value_at(i);

//             std::cout << "k: " << key
//                 << " -> (" << value.first << ", "
//                 << value.second << ")\n";

//           }
//         });
//   }
//   Kokkos::finalize();

//   return 0;
// }