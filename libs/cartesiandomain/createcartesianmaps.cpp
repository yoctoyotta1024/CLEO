/*
 * ----- CLEO -----
 * File: createcartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
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

void set_1Dmodel_maps(const GbxBoundsFromBinary &gfb,
                      CartesianMaps &gbxmaps);

void set_1Dmodel_areas_vols(const GbxBoundsFromBinary &gfb,
                            CartesianMaps &gbxmaps);

bool at_domainboundary(const unsigned int idx,
                       const unsigned int increment,
                       const unsigned int ndim);

Kokkos::pair<unsigned int, unsigned int>
finitedomain_nghbours(const unsigned int idx,
                      const unsigned int increment,
                      const unsigned int ndim);

Kokkos::pair<unsigned int, unsigned int>
cartesian_znghbrs(const unsigned int idx,
                  const std::vector<size_t> &ndims);

Kokkos::pair<double, double> nullbounds()
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
  const GbxBoundsFromBinary gfb(nspacedims, grid_filename);

  CartesianMaps gbxmaps(gfb.get_ngbxs());

  set_maps_ndims(gfb.ndims, gbxmaps);

  if (nspacedims == 0)
  {
    set_0Dmodel_maps(gfb, gbxmaps);
    set_0Dmodel_areas_vols(gfb, gbxmaps);
  }

  // else if (nspacedims == 1)
  // {
  //   set_1Dmodel_maps(gfb);
  //   set_1Dmodel_areas_vols(gfb, gbxmaps);
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
  viewd_ndims d_ndims("ndims");                       // view for ndims (on device)
  auto h_ndims = Kokkos::create_mirror_view(d_ndims); // mirror ndims in case view is on device memory

  for (unsigned int m(0); m < 3; ++m)
  {
    h_ndims(m) = ndims.at(m);
  }
  Kokkos::deep_copy(d_ndims, h_ndims);

  gbxmaps.set_ndims(d_ndims);
}

void set_0Dmodel_maps(const GbxBoundsFromBinary &gfb,
                      CartesianMaps &gbxmaps)
/* gives all coord[X]bounds maps to 1 key with null
values (max/min numerical limits). Sets periodic
boundary conditions in all directions meaning neighbour
of single gridbox with gbxidx=0 is itself */
{
  gbxmaps.to_coord3bounds.insert(0, nullbounds());
  gbxmaps.to_coord1bounds.insert(0, nullbounds());
  gbxmaps.to_coord2bounds.insert(0, nullbounds());

  gbxmaps.to_back_coord3nghbr.insert(0, nullnghbr(0));
  gbxmaps.to_forward_coord3nghbr.insert(0, nullnghbr(0));
  gbxmaps.to_back_coord1nghbr.insert(0, nullnghbr(0));
  gbxmaps.to_forward_coord1nghbr.insert(0, nullnghbr(0));
  gbxmaps.to_back_coord2nghbr.insert(0, nullnghbr(0));
  gbxmaps.to_forward_coord2nghbr.insert(0, nullnghbr(0));
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

void set_1Dmodel_maps(const GbxBoundsFromBinary &gfb,
                      CartesianMaps &gbxmaps)
/* Gives all coord[X]bounds maps for X = x or y
null values (max/min numerical limits) for all gridboxes and
periodic boundary conditions (meaning neighbour
of gridbox in x or y is itself). coord3bounds, ie. z direction
maps are set using vecors from gfb. It is assumed that for a
gridbox with it's index at position 'p' in the gfb.gbxidxs
vector, the [zmin, zmax] coords of that gridbox are at
[pos, pos+1] in the gfb.gbxidxs vector, where pos = p*6 */
{
  for (auto idx : gfb.gbxidxs)
  {
    gbxmaps.to_coord3bounds.insert(idx, gfb.get_gbxbounds(idx));
    gbxmaps.to_coord1bounds.insert(idx, nullbounds());
    gbxmaps.to_coord2bounds.insert(idx, nullbounds());

    // gbxmaps.to_back_coord3nghbr.insert(0, nullnghbr(0));
    // gbxmaps.to_forward_coord3nghbr.insert(0, nullnghbr(0));
    gbxmaps.to_back_coord1nghbr.insert(idx, nullnghbr(0));
    gbxmaps.to_forward_coord1nghbr.insert(idx, nullnghbr(0));
    gbxmaps.to_back_coord2nghbr.insert(idx, nullnghbr(0));
    gbxmaps.to_forward_coord2nghbr.insert(idx, nullnghbr(0));
  }
}

void set_1Dmodel_areas_vols(const GbxBoundsFromBinary &gfb,
                            CartesianMaps &gbxmaps)
/* sets (finite) dimensionless horizontal area and volume
of single gridbox in 0-D model (ie. entire domain) */
{
  // const double domainarea = gfb.gbxarea_fromgridfile(0);
  // const double domainvol = gfb.gbxvol_fromgridfile(0);

  // gbxmaps.set_gbxarea(domainarea);
  // gbxmaps.set_gbxvolume(domainvol);
}

bool at_domainboundary(const unsigned int idx,
                       const unsigned int increment,
                       const unsigned int ndim)
/* returns true if idx for gridbox is at a domain boundary, given
neighbouring indexes are +- increment from idx and the number of
gridboxes making up the domain in that direction (ndim) */
{
  return (idx / increment) % ndim == 0;
}

Kokkos::pair<unsigned int, unsigned int>
finitedomain_nghbours(const unsigned int idx,
                      const unsigned int increment,
                      const unsigned int ndim)
/* returns {forward, backward} gridbox neighbours with
treatment of neighbours as if bounds of domain are finite.
This means that no neighbour exists above highest / below lowest
gridboxes in a given direction. Non-existent neighbours are
defined with gbxindex = max unsigned int, meaning in a given
direction the gbxindex of the backwards / forwards neighbour
of a gridbox at the edge of the domain is a max unsigned int */
{
  unsigned int forward(idx + increment);
  unsigned int backward(idx - increment);

  if (at_domainboundary(idx, increment, ndim)) // at lower edge of domain
  {
    backward = LIMITVALUES::uintmax;
  }

  if (at_domainboundary(forward, increment, ndim)) // at upper edge of domain
  {
    forward = LIMITVALUES::uintmax;
  }

  return {forward, backward};
}

Kokkos::pair<unsigned int, unsigned int>
cartesian_znghbrs(const unsigned int idx,
                  const std::vector<size_t> &ndims)
/* returns pair of gbx indexes for {backwards, forwards}
neighbour of a gridbox with index 'idx' in cartesian domain.
Treatment of neighbours for gridboxes at the edges of the
domain is either finite (null neighbour) or periodic
(cyclic neighbour) */
{
  return finitedomain_nghbours(idx, 1, ndims.at(0));
  // return periodicdomain_nghbours(idx, 1, ndims.at(0));
}