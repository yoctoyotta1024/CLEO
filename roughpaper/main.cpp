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
using kkpair = Kokkos::pair<size_t, size_t>;


kkpair set_refs(const unsigned int ii, viewd_constsupers supers)
{
  const kkpair new_refs = {0,0};
  return new_refs;
}


int main(int argc, char *argv[])
{
  const size_t nsupers(20);
  const size_t ngbxs(1);

  Kokkos::initialize(argc, argv);
  {
    viewd_supers supers("supers", nsupers);
    auto h_supers = Kokkos::create_mirror_view(supers);
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      const unsigned int ii(kk);
      h_supers(kk) = Superdrop();
      h_supers(kk).set_sdgbxindex(ii);
      std::cout << "ii: " << h_supers(kk).get_sdgbxindex() << "\n";
    }
    Kokkos::deep_copy(supers, h_supers);

    Kokkos::View<kkpair *> viewd_refs("gbxs", ngbxs); 
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      kkpair refs{0, 0};
      refs = set_refs(ii, supers);
      std::cout << "refs: " << refs.first << ", " << refs.second << "\n";
      
      viewd_refs(ii) = refs;
    }

    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      std::cout << "refs: " << viewd_refs(ii).first
                << ", " << viewd_refs(ii).second << "\n";
      const auto subview = Kokkos::subview(supers, viewd_refs(ii));
      const size_t n(subview.extent(0));
      std::cout << "---- gbx: ii = " << ii << " -----\n";

      for (size_t kk(0); kk < n; ++kk)
      {
        std::cout << "ii: " << subview(kk).get_sdgbxindex() << "\n"; 
      }
      std::cout << "n = " << n << "\n"; 
    }
  }
  Kokkos::finalize();
}