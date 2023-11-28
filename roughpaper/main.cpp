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
private:
  unsigned int sdgbxindex;
public:

KOKKOS_INLINE_FUNCTION
void set_sdgbxindeX(const unsigned int ii){sdgbxindex = ii}

KOKKOS_INLINE_FUNCTION
    unsigned int get_sdgbxindeX() const
{
  return sdgbxindex;
}
};

using viewd_supers = Kokkos::View<Superdrop *>;
using ExecSpace = Kokkos::DefaultExecutionSpace;

int main(int argc, char *argv[])
{
  const size_t nsupers(20);

  Kokkos::initialize(argc, argv);
  {
    viewd_supers supers("supers", nsupers);
    auto h_supers = Kokkos::create_mirror_view(supers); // mirror of supers in case view is on device memory
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      const unsigned int ii(kk/2);
      h_supers(kk) = Superdrop{ii};
      std::cout << "ii: " << h_supers(kk).get_sdgbxindex() << "\n";
    }
    Kokkos::deep_copy(supers, h_supers);
     
  }
  Kokkos::finalize();
}