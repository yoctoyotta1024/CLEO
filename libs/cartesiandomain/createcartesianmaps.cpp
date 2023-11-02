/*
 * ----- CLEO -----
 * File: createcartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 2nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for creating a cartesian maps struct
 * from a GbxBoundsFromBinary struct containing
 * vectors of gridbox's indexes and their
 * coordinate (upper and lower) boundaries
 */

#include "./createcartesianmaps.hpp"

void set_maps_ndims(const std::vector<size_t> &ndims,
                    CartesianMaps &gbxmaps);

void set_0Dmodel_maps(const GbxBoundsFromBinary &gfb,
                      CartesianMaps &gbxmaps);

void set_0Dmodel_areas_vols(const GbxBoundsFromBinary &gfb,
                            CartesianMaps &gbxmaps);

std::pair<double, double> nullbounds()
/* bounds for CartesianMaps of gridboxes along
directions of model not used e.g. in 1-D model,
these are bounds of gridboxes in coord1 and
coord2 directions */
{
  return {LIMITVALUES::llim, LIMITVALUES::ulim};
}

unsigned int nullnghbr(const unsigned int idx)
/* neighbours for CartesianMaps of gridboxes along
directions of model not used. Boundaries are
'periodic' BCs in non-existent dimensions 
e.g. in 2-D model, neighbour in coord2 direction
of gridbox is itself */
{
  return idx;
}

CartesianMaps create_cartesian_maps(const unsigned int nspacedims,
                                    std::string_view grid_filename)
/* creates cartesian maps instance using gridbox bounds read from 
gridfile. In a non-3D case, boundaries for unused dimensions may be
the min/max possible (numerical limits), however the area and volume
of each gridbox remains finite. E.g. In the 0-D case, the maps have 1
{key, value} for gridbox 0 which are numerical limits, whilst the
volume function returns a value determined from the gridfile input */
{
  CartesianMaps gbxmaps;
  const GbxBoundsFromBinary gfb(nspacedims, grid_filename);

  set_maps_ndims(gfb.ndims, gbxmaps);
  
  if (nspacedims == 0)
  {
    set_0Dmodel_maps(gfb, gbxmaps);
    set_0Dmodel_areas_vols(gfb, gbxmaps);
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

  // else
  // {
  //   throw std::invalid_argument("nspacedims > 3 is invalid ");
  // }

  // check_ngridboxes();

  return CartesianMaps();
}

void set_maps_ndims(const std::vector<size_t> &ndims,
                    CartesianMaps &gbxmaps)
/* copys ndims  to gbxmaps' ndims to set number of
dimensions (ie. number of gridboxes) in
[coord3, coord1, coord2] directions */
{
  viewd_ndims d_ndims("ndims"); // view for ndims (on device)
  auto h_ndims = Kokkos::create_mirror_view(d_ndims); // mirror ndims in case view is on device memory

  for (unsigned int m(0); m < 3; ++m)
  {
    h_ndims(m)= ndims.at(m);
  }
  Kokkos::deep_copy(d_ndims, h_ndims);

  gbxmaps.set_ndims(d_ndims);
}

void set_0Dmodel_maps(const GbxBoundsFromBinary &gfb,
                      CartesianMaps &gbxmaps)
/* set all coord[X]bounds maps to null values
(max/min numerical limits). Sets periodic boundary
conditions in all directions meaning neighbour of
single gridbox with gbxidx=0 is itself */
{ 
  kokkos_pairmap 
  gbxmaps.set_boundsmaps_via_copy(nullbounds(),
                                  nullbounds(),
                                  nullbounds())
  
  gbxmaps.set_nghbrsmaps_via_copy(3, nullnghbr(0), nullnghbr(0));
  gbxmaps.set_nghbrsmaps_via_copy(1, nullnghbr(0), nullnghbr(0));
  gbxmaps.set_nghbrsmaps_via_copy(2, nullnghbr(0), nullnghbr(0));
}

void set_0Dmodel_areas_vols(const GbxBoundsFromBinary &gfb,
                            CartesianMaps &gbxmaps)
/* sets (finite) dimensionless horizontal area and volume
of single gridbox in 0-D model (ie. entire domain) */
{
  const double domainarea = gfb.gbxarea_fromgridfile(0);
  const double domainvol = gfb.gbxvol_fromgridfile(0);

  gbxmaps.set_gbxarea(domainarea);
  gbxmaps.set_gbxvolume(domainvol);
}