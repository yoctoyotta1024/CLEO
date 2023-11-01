/*
 * ----- CLEO -----
 * File: createcartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 1st November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for creating a cartesian maps struct
 * from a GridBoxBoundaries struct containing gridbox's
 * indexes and their coordinate (upper and lower) boundaries
 */

#include "./createmaps_frombinary.hpp"

void set_maps_ndims(const std::vector<size_t> &ndims,
                    CartesianMaps &gbxmaps);

void set_0Dmodel_maps();

void set_0Dmodel_gbxvolumes();

CartesianMaps create_cartesian_maps(const unsigned int nspacedims,
                                    std::string_view grid_filename)
{
  CartesianMaps gbxmaps;
  const GbxBoundsFromBinary gfb(nspacedims, grid_filename);

  set_maps_ndims(gfb.ndims, gbxmaps);

  if (nspacedims == 0)
  {
    set_0Dmodel_maps();
    set_0Dmodel_gbxvolumes();
  }

  // else if (nspacedims == 1)
  // {
  //   set_1Dmodel_maps(gfb);
  // }

  // else if (nspacedims == 2)
  // {
  //   set_2Dmodel_maps(gfb);
  // }

  // else if (nspacedims == 3)
  // {
  //   set_3Dmodel_maps(gfb);
  // }

  else
  {
    throw std::invalid_argument("nspacedims > 3 is invalid ");
  }

  check_ngridboxes();

  return CartesianMaps();
}

void set_maps_ndims(const std::vector<size_t> &ndims,
                    CartesianMaps &gbxmaps)
/* copys ndims  to gbxmaps' ndims to set number of
dimensions (ie. number of gridboxes) in
[coord3, coord1, coord2] directions */
{
  auto h_ndims = Kokkos::create_mirror_view(gbxmaps.ndims); // mirror ndims in case view is on device memory

  for (unsigned int m(0); m < 3; ++m)
  {
    h_ndims(m) = ndims.at(m);
  }

  Kokkos::deep_copy(gbxmaps.ndims, h_ndims);
}

void set_0Dmodel_maps(const double domainarea,
                      const double domainvol)
/* set idx2bounds_[i] maps to numeical limits. Set volume
 map using coords read from gridfile */
{
  idx2bounds_z[0] = numeric_limit_bounds();
  idx2bounds_x[0] = numeric_limit_bounds();
  idx2bounds_y[0] = numeric_limit_bounds();
  
  idx2nghbour_z[0] = {0, 0}; // 'periodic' BCs in non-existent dimensions 
  idx2nghbour_x[0] = {0, 0};
  idx2nghbour_y[0] = {0, 0};
}

void set_0Dmodel_gbxvolumes()
{
  const double domainarea = get_0Ddomainarea_from_gridfile(gfb);
  const double domainvol = get_0Ddomainvol_from_gridfile(gfb);

  idx2area[0] = domainarea; // dimensionless horizontal area of 0D model
  idx2vol[0] = domainvol; // dimensionless volume of 0D model
}