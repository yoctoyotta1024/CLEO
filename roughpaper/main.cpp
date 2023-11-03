/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
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

#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>

void set_0Dmodel_maps()
{
  std::cout << "0-D\n";
}

void set_1Dmodel_maps()
{
  std::cout << "1-D\n";
}

int main(int argc, char *argv[])
{
  const unsigned int nspacedims(3);

  if (nspacedims == 0)
  {
    set_0Dmodel_maps();
  }

  else if (nspacedims == 1)
  {
    set_1Dmodel_maps();
  }

  switch (nspacedims)
  {
  case 0:
    set_0Dmodel_maps();
    break;
  case 1:
    set_1Dmodel_maps();
    break;
  default:
    std::cout << "error\n";
  }
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