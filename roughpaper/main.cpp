/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Wednesday 1st November 2023
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
 * rough paper for checking small things
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <limits>

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <Kokkos_Random.hpp>

using viewd_supers = Kokkos::View<int *>;

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

int main(int argc, char *argv[])
{
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using GenRandomPool = Kokkos::Random_XorShift64_Pool<ExecSpace>; // type for pool of thread safe random number generators

  std::vector<int> i_supers = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  size_t nsupers(i_supers.size());

  Kokkos::initialize(argc, argv);
  {
    viewd_supers supers("supers", nsupers);

    auto h_supers = Kokkos::create_mirror_view(supers); // mirror of supers in case view is on device memory
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      h_supers(kk) = i_supers.at(kk);
    }
    Kokkos::deep_copy(supers, h_supers);

    GenRandomPool genpool(std::random_device{}());
    {
      URBG<ExecSpace> urbg{genpool.get_state()}; // thread safe random number generator

      std::cout << " \n --- b4 ---\n ";
      for (size_t kk(0); kk < nsupers; ++kk)
      {
        std::cout << supers(kk) << ", ";
      }
      std::cout << " \n --- --- ---\n ";

      supers = shuffle_supers(supers, urbg);

      std::cout << " \n --- l8r ---\n ";
      for (size_t kk(0); kk < nsupers; ++kk)
      {
        std::cout << supers(kk) << ", ";
      }
      std::cout << " \n --- --- ---\n ";

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