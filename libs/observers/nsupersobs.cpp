/*
 * ----- CLEO -----
 * File: nsupersobs.cpp
 * Project: observers
 * Created Date: Thursday 30th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 30th November 2023 
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality to output nsupers
 * (per gridbox or total in domain)
 * to array in a zarr file system storage
 */

#include "./nsupersobs.hpp"

size_t calc_nrainsupers(const SupersInGbx &supersingbx)
/* returns count of number of "raindrop-like" superdrops
for each gridbox. "raindrop-like" means radius > rlim.
  * WARNING! * When using OpenMP (supers in Host Space)
 and there are only a few superdroplets in supers,
 calc_nrainsupers is much slower then calc_nrainsupers_serial
 (probably because opening threads is more costly than the
 time saved in a parallel calculation over few elements) */
{
  constexpr double rlim(40e-6 / dlc::R0); // dimless minimum radius of raindrop
  
  const auto h_supers = supersingbx.hostcopy();

  size_t nrainsupers(0);
  for (size_t kk(0); kk < h_supers.extent(0); ++kk)
  {
    
  }

  // Kokkos::parallel_reduce(
  //     "calc_massmoments",
  //     Kokkos::RangePolicy<ExecSpace>(0, nsupers),
  //     KOKKOS_LAMBDA(const size_t kk, size_t &nrs) {
  //       const double xi = (double)(d_supers(kk).get_xi()); // cast multiplicity from unsigned int to double

  //       const auto radius = d_supers(kk).get_radius();
  //       if (radius >= rlim)
  //       {
  //         ++nrainsupers;
  //       }
  //     },
  //     nrainsupers); // {0th, 1st, 2nd} mass moments

  return nrainsupers;
}

size_t calc_nrainsupers_serial(const SupersInGbx &supersingbx)
/* deep copy if necessary (if superdrops are on device not
  host memory), then returns count of number of "raindrop-like" 
  superdrops for each gridbox. "raindrop-like" means radius > rlim */
{
  constexpr double rlim(40e-6 / dlc::R0); // dimless minimum radius of raindrop
  
  const auto h_supers = supersingbx.hostcopy();

  size_t nrainsupers(0);
  for (size_t kk(0); kk < h_supers.extent(0); ++kk)
  {
    const auto radius = h_supers(kk).get_radius();
    if (radius >= rlim)
    {
      ++nrainsupers;
    }
  }

  return nrainsupers;
}