/*
 * ----- CLEO -----
 * File: massmomentsobs.cpp
 * Project: observers
 * Created Date: Sunday 22nd October 2023
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
 * functionality to calculate and output mass moments 
 * of droplet distribution to array in a zarr file
 * system storage
 */

#include "./massmomentsobs.hpp"

void DoMassMomentsObs::
    massmoments_to_storage(const mirrorh_constsupers h_supers) const
/* calculated 0th, 1st and 2nd moment of the (real) droplet mass
distribution and then writes them to zarr storage. (i.e.
0th, 3rd and 6th moment of the droplet radius distribution) */
{
  std::array<double, 3> moms({0.0, 0.0, 0.0}); // 0th, 1st and 2nd mass moments
  for (size_t kk(0); kk < h_supers.extent(0); ++kk)
  {
    const double xi = (double)(h_supers(kk).get_xi()); // cast multiplicity from unsigned int to double
    const double mass(h_supers(kk).mass());
    moms.at(0) += xi;
    moms.at(1) += xi * mass;
    moms.at(2) += xi * mass * mass;
  }

  zarr->values_to_storage(moms);
}

void DoRainMassMomentsObs::
    rainmassmoments_to_storage(const mirrorh_constsupers h_supers) const
/* calculated 0th, 1st and 2nd moment of the (real) droplet mass
distribution and then writes them to zarr storage. (i.e.
0th, 3rd and 6th moment of the droplet radius distribution) */
{
  constexpr double rlim(40e-6 / dlc::R0); // dimless minimum radius of raindrop

  std::array<double, 3> moms({0.0, 0.0, 0.0}); // 0th, 1st and 2nd mass moments
  for (size_t kk(0); kk < h_supers.extent(0); ++kk)
  {
    if (h_supers(kk).get_radius() >= rlim)
    {
      const double xi = (double)(h_supers(kk).get_xi()); // cast multiplicity from unsigned int to double
      const double mass(h_supers(kk).mass());
    
      moms.at(0) += xi;
      moms.at(1) += xi * mass;
      moms.at(2) += xi * mass * mass;
    }
  }

  zarr->values_to_storage(moms);
}