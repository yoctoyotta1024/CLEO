/*
 * ----- CLEO -----
 * File: massmomentsobs.cpp
 * Project: observers
 * Created Date: Sunday 22nd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 23rd October 2023
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

std::array<double, 3>
calc_massmoments(const subviewd_constsupers d_supers)
/* calculated 0th, 1st and 2nd moment of the (real)
droplet mass distribution, i.e. 0th, 3rd and 6th
moment of the droplet radius distribution.
Kokkos::parallel_reduce([...]) is equivalent in serial to:
for (size_t kk(0); kk < d_supers.extent(0); ++kk){[...]},
see calc_massmoments_serial */
{
  std::array<double, 3> moms({0.0, 0.0, 0.0}); // {0th, 1st, 2nd} mass moments

  const size_t nsupers(d_supers.extent(0));
  Kokkos::parallel_reduce(
      "massmoments_to_storage",
      Kokkos::RangePolicy<ExecSpace>(0, nsupers),
      KOKKOS_LAMBDA(const size_t kk, double &m0, double &m1, double &m2) {
        
        const double xi = (double)(d_supers(kk).get_xi()); // cast multiplicity from unsigned int to double
        const double mass(d_supers(kk).mass()); 
        m0 += xi;
        m1 += xi * mass;
        m2 += xi * mass * mass;
      },
      moms.at(0), moms.at(1), moms.at(2)); // {0th, 1st, 2nd} mass moments

  return moms;
}

std::array<double, 3>
calc_rainmassmoments(const subviewd_constsupers d_supers)
/* calculated 0th, 1st and 2nd moment of the
(real) raindroplet mass distribution, i.e. 0th, 3rd and 6th
moment of the droplet radius disttribution. Raindrops are 
all droplets with r >= rlim = 40 microns.
Kokkos::parallel_reduce([...]) is equivalent in serial to:
for (size_t kk(0); kk < d_supers.extent(0); ++kk){[...]},
see calc_rainmassmoments_serial */
{
  constexpr double rlim(40e-6 / dlc::R0); // dimless minimum radius of raindrop

  std::array<double, 3> moms({0.0, 0.0, 0.0}); // {0th, 1st, 2nd} rain mass moments

  const size_t nsupers(d_supers.extent(0));
  Kokkos::parallel_reduce(
      "rainmassmoments_to_storage",
      Kokkos::RangePolicy<ExecSpace>(0, nsupers),
      KOKKOS_LAMBDA(const size_t kk, double &m0, double &m1, double &m2) {
        if (h_supers(kk).get_radius() >= rlim)
        {
          const double xi = (double)(d_supers(kk).get_xi()); // cast multiplicity from unsigned int to double
          const double mass(d_supers(kk).mass());
          m0 += xi;
          m1 += xi * mass;
          m2 += xi * mass * mass;
        }
      },
      moms.at(0), moms.at(1), moms.at(2)); // {0th, 1st, 2nd} rain mass moments

  return moms;
}

std::array<double, 3>
calc_massmoments_serial(const subviewd_constsupers d_supers)
/* calculated 0th, 1st and 2nd moment of the (real)
droplet mass distribution, i.e. 0th, 3rd and 6th
moment of the droplet radius distribution */
{
  std::array<double, 3> moms({0.0, 0.0, 0.0}); // {0th, 1st, 2nd} mass moments

  auto h_supers = Kokkos::create_mirror_view(d_supers);
  Kokkos::deep_copy(h_supers, d_supers);
  for (size_t kk(0); kk < h_supers.extent(0); ++kk)
  {
    const double xi = (double)(h_supers(kk).get_xi()); // cast multiplicity from unsigned int to double
    const double mass(h_supers(kk).mass());
    moms.at(0) += xi;
    moms.at(1) += xi * mass;
    moms.at(2) += xi * mass * mass;
  }

  return moms;
}

std::array<double, 3>
calc_rainmassmoments_serial(const subviewd_constsupers d_supers)
/* calculated 0th, 1st and 2nd moment of the (real)
raindroplet mass distribution, i.e. 0th, 3rd and 6th
moment of the droplet radius disttribution. Raindrops
are all droplets with r >= rlim = 40 microns */
{
  std::array<double, 3> moms({0.0, 0.0, 0.0}); // {0th, 1st, 2nd} mass moments

  auto h_supers = Kokkos::create_mirror_view(d_supers);
  Kokkos::deep_copy(h_supers, d_supers);
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

  return moms;
}

void DoMassMomentsObs::
    massmoments_to_storage(const subviewd_constsupers d_supers) const
/* calculated 0th, 1st and 2nd moment of the (real) droplet mass
distribution and then writes them to zarr storage. (I.e.
0th, 3rd and 6th moment of the droplet radius distribution).
Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t kk(0); kk < d_supers.extent(0); ++kk){[...]} */
{
  const auto moms = calc_massmoments(d_supers);

  zarr->values_to_storage(moms); // {0th, 1st, 2nd} mass moments
}

void DoRainMassMomentsObs::
    rainmassmoments_to_storage(const mirrorh_constsupers h_supers) const
/* calculated 0th, 1st and 2nd moment of the (real) droplet mass
distribution and then writes them to zarr storage. (I.e.
0th, 3rd and 6th moment of the droplet radius distribution) */
{
  const auto moms = calc_rain_massmoments(d_supers);

  zarr->values_to_storage(moms); // {0th, 1st, 2nd} rain mass moments
}