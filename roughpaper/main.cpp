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

viewd_supers init_supers(const size_t nsupers,
                         const size_t ngbxs)
{
    viewd_supers supers("supers", nsupers);
    auto h_supers = Kokkos::create_mirror_view(supers);
    for (size_t kk(0); kk < nsupers; ++kk)
    {
      const unsigned int ii(1);
      h_supers(kk) = Superdrop();
      h_supers(kk).set_sdgbxindex(ii);
      std::cout << "ii: " << h_supers(kk).get_sdgbxindex() << "\n";
    }
    Kokkos::deep_copy(supers, h_supers);

  return supers;
}

unsigned int make_obs_gbx(unsigned int nobs,
                          unsigned int nobs_gbxs,
                          const size_t ngbxs)
{

  ++nobs_gbxs;

  std::cout << "nobs: " << nobs
            << " & " << nobs_gbxs << " -> "
            << nobs_gbxs / ngbxs << "\n";

  nobs = nobs_gbxs / ngbxs;

  return nobs;
}

unsigned int make_obs_gbxs(unsigned int nobs)
{
  return ++nobs;
}

int main(int argc, char *argv[])
{
  const size_t nsupers(10);
  const size_t ngbxs(13);
  const size_t nobs(1);

  Kokkos::initialize(argc, argv);
  {
    auto supers = init_supers(nsupers, ngbxs);

    unsigned int nobs_1(0);
    unsigned int nobs_2(0);
    unsigned int nobs_gbxs(0);

    for (size_t jj(0); jj < nobs; ++jj)
    {

      for (size_t ii(0); ii < ngbxs; ++ii)
      {
        nobs_1 = make_obs_gbx(nobs_1, nobs_gbxs, ngbxs);
        ++nobs_gbxs;
      }

      nobs_2 = make_obs_gbxs(nobs_2);
    }

    std::cout << "final nobs: " << nobs_1
              << " == " << nobs_2 << "\n";
  }
  Kokkos::finalize();
}