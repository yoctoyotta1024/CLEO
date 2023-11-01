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

CartesianMaps create_cartesian_maps(const unsigned int nspacedims,
                                    std::string_view grid_filename)
{
  CartesianMaps gbxmaps;
  const GbxBoundsFromBinary gfb(nspacedims, grid_filename);

  set_maps_ndims(gfb.ndims, gbxmaps);

  // if (nspacedims == 0)
  // {
  //   const double domainarea = get_0Ddomainarea_from_gridfile(gfb);
  //   const double domainvol = get_0Ddomainvol_from_gridfile(gfb);
  //   set_0Dmodel_maps(domainarea, domainvol);
  // }

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

  // else
  // {
  //   throw std::invalid_argument("nspacedims > 3 is invalid ");
  // }

  // check_ngridboxes();

  // return CartesianMaps();
}

void set_maps_ndims(const std::vector<size_t> &ndims,
                    CartesianMaps &gbxmaps)
/* sets dimensions (ie. number of gridboxes)
in [coord3, coord1, coord2] directions */
{
  auto h_ndims = Kokkos::create_mirror_view(gbxmaps.ndims); // mirror ndims in case view is on device memory

  for (unsigned int m(0); m < 3; ++m)
  {
    h_ndims(m) = ndims.at(m);
  }

  Kokkos::deep_copy(gbxmaps.ndims, h_ndims);
}