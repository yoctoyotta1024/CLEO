/*
 * ----- CLEO -----
 * File: main.cpp
 * Project: roughpaper
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 18th December 2023
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

// ---------------------------------------------------------------- //

struct Superdrop
{
private:
  unsigned int sdgbxindex;
public:

KOKKOS_INLINE_FUNCTION
void set_sdgbxindex(const unsigned int ii){sdgbxindex = ii;}

KOKKOS_INLINE_FUNCTION
    unsigned int get_sdgbxindex() const
{
  return sdgbxindex;
}
};

using viewd_supers = Kokkos::View<Superdrop *>;
using viewd_constsupers = Kokkos::View<const Superdrop *>; // view in device memory of const superdroplets
using ExecSpace = Kokkos::DefaultExecutionSpace;

viewd_supers init_supers(const size_t nsupers,
                         const size_t ngbxs)
{
    viewd_supers supers("supers", nsupers);
    auto h_supers = Kokkos::create_mirror_view(supers);
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      const auto ii = unsigned int(1);
      h_supers(kk) = Superdrop();
      h_supers(kk).set_sdgbxindex(ii);
      std::cout << "ii: " << h_supers(kk).get_sdgbxindex() << "\n";
    }
    Kokkos::deep_copy(supers, h_supers);

  return supers;
}

// ---------------------------------------------------------------- //

int main(int argc, char *argv[])
{
  const size_t nsupers(10);
  const size_t ngbxs(13);

  Kokkos::initialize(argc, argv);
  {
    auto supers = init_supers(nsupers, ngbxs);

    for (size_t ii(0); ii < ngbxs; ++ii)
    {
    }
  }
  
  Kokkos::finalize();
}