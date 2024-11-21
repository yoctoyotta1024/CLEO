/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: createcartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 1st May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for creating a cartesian maps struct
 * from a GbxBoundsFromBinary struct containing
 * vectors of gridbox's indexes and their
 * coordinate (upper and lower) boundaries
 */

#include "cartesiandomain/createcartesianmaps.hpp"

#include "mpi.h"

void check_ngridboxes_matches_ndims(const CartesianMaps &gbxmaps, const size_t ngbxs);

void check_ngridboxes_matches_maps(const CartesianMaps &gbxmaps, const size_t ngbxs);

void set_maps_ndims(const std::vector<size_t> &ndims, CartesianMaps &gbxmaps);

void set_model_areas_vols(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps);

/* set gbxindex to out of bounds value */
void set_outofbounds(CartesianMaps &gbxmaps);

void set_0Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps);

void set_1Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps);

void set_2Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps);

void set_3Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps);

/* bounds for CartesianMaps of gridboxes along
directions of model not used e.g. in 1-D model,
these are bounds of gridboxes in coord1 and
coord2 directions */
Kokkos::pair<double, double> nullbounds() { return {LIMITVALUES::llim, LIMITVALUES::ulim}; }

/* {back, forward} neighbours for CartesianMaps of
gridboxes along directions of model not used. Boundaries
are 'periodic' BCs in non-existent dimensions
e.g. in 2-D model, neighbour in coord2 direction
of gridbox is itself */
Kokkos::pair<unsigned int, unsigned int> nullnghbrs(const unsigned int idx) { return {idx, idx}; }

/* creates cartesian maps instance using gridbox bounds read from
gridfile for a 0-D, 1-D, 2-D or 3-D model with periodic or finite
boundary conditions. In a non-3D case, boundaries and neighbours
maps for unused dimensions are 'null' (ie. return numerical limits),
however the area and volume of each gridbox remains finite.
E.g. In the 0-D case, the bounds maps all have 1 {key, value} where
key=gbxidx=0 and value = {max, min} numerical limits, meanwhile volume
function returns a value determined from the gridfile 'grid_filename' */
CartesianMaps create_cartesian_maps(const size_t ngbxs, const unsigned int nspacedims,
                                    const std::filesystem::path grid_filename) {
  std::cout << "\n--- create cartesian gridbox maps ---\n";

  const GbxBoundsFromBinary gfb(ngbxs, nspacedims, grid_filename);

  CartesianMaps gbxmaps(gfb.get_ngbxs());

  gbxmaps.create_decomposition(
      gfb.ndims, gfb.get_coord3gbxbounds(0).second - gfb.get_coord3gbxbounds(0).first,
      gfb.get_coord1gbxbounds(0).second - gfb.get_coord1gbxbounds(0).first,
      gfb.get_coord2gbxbounds(0).second - gfb.get_coord2gbxbounds(0).first);

  set_maps_ndims(gfb.ndims, gbxmaps);
  set_outofbounds(gbxmaps);

  switch (nspacedims) {
    case 0:
      set_model_areas_vols(gfb, gbxmaps);
      set_0Dmodel_maps(gfb, gbxmaps);
      break;
    case 1:
      set_model_areas_vols(gfb, gbxmaps);
      set_1Dmodel_maps(gfb, gbxmaps);
      break;
    case 2:
      set_model_areas_vols(gfb, gbxmaps);
      set_2Dmodel_maps(gfb, gbxmaps);
      break;
    case 3:
      set_3Dmodel_maps(gfb, gbxmaps);
      break;
    default:
      throw std::invalid_argument("nspacedims > 3 is invalid ");
  }

  check_ngridboxes_matches_ndims(gbxmaps, gbxmaps.get_total_global_ngridboxes());
  check_ngridboxes_matches_maps(gbxmaps, gbxmaps.get_local_ngridboxes_hostcopy());

  std::cout << "--- create cartesian gridbox maps: success ---\n";

  return gbxmaps;
}

/* checks number of gridboxes according to
maps matches with expected value from gfb */
void check_ngridboxes_matches_ndims(const CartesianMaps &gbxmaps, const size_t ngbxs) {
  const auto h_ndims = gbxmaps.get_ndims_hostcopy();
  const size_t ngbxs_from_ndims(h_ndims(0) * h_ndims(1) * h_ndims(2));

  if (ngbxs_from_ndims != ngbxs) {
    throw std::invalid_argument(
        "ndims from gridbox maps inconsistent "
        " with number of gridboxes");
  }
}

void check_ngridboxes_matches_maps(const CartesianMaps &gbxmaps, const size_t ngbxs) {
  const size_t ngbxs_from_maps(gbxmaps.maps_size());
  if (ngbxs_from_maps != ngbxs + 1) {
    throw std::invalid_argument(
        "ngbxs from gridbox maps inconsistent "
        " with number of gridboxes");
  }
}

/* copys ndims  to gbxmaps' ndims to set number of
dimensions (ie. number of gridboxes) in
[coord3, coord1, coord2] directions */
void set_maps_ndims(const std::vector<size_t> &i_ndims, CartesianMaps &gbxmaps) {
  auto h_ndims = Kokkos::create_mirror_view(
      gbxmaps.get_ndims());  // mirror ndims in case view is on device memory

  for (unsigned int m(0); m < 3; ++m) {
    h_ndims(m) = i_ndims.at(m);
  }

  gbxmaps.set_ndims_via_copy(h_ndims);
}

/* sets (finite) dimensionless horizontal area and
volume using area and volume from gfb for gbxidx=0 */
void set_model_areas_vols(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps) {
  for (auto idx : gfb.gbxidxs) {
    gbxmaps.insert_gbxarea(idx, gfb.gbxarea(idx));
    gbxmaps.insert_gbxvolume(idx, gfb.gbxvol(idx));
  }
}

/* sets value for coordinate bounds for case
when outofbounds gbxidx searches map */
void set_outofbounds(CartesianMaps &gbxmaps) {
  const auto idx = (unsigned int)LIMITVALUES::oob_gbxindex;
  gbxmaps.insert_coord3bounds(idx, nullbounds());
  gbxmaps.insert_coord1bounds(idx, nullbounds());
  gbxmaps.insert_coord2bounds(idx, nullbounds());

  gbxmaps.insert_coord3nghbrs(idx, nullnghbrs(idx));
  gbxmaps.insert_coord1nghbrs(idx, nullnghbrs(idx));
  gbxmaps.insert_coord2nghbrs(idx, nullnghbrs(idx));

  gbxmaps.insert_gbxarea(idx, 0.0);
  gbxmaps.insert_gbxvolume(idx, 0.0);
}

/* gives all coord[X]bounds maps to 1 key with null
values (max/min numerical limits). Also sets null
neighbours maps (meaning periodic boundary conditions
in all directions where neighbour of single gridbox
with gbxidx=0 is itself) */
void set_0Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps) {
  gbxmaps.insert_coord3bounds(0, nullbounds());
  gbxmaps.insert_coord1bounds(0, nullbounds());
  gbxmaps.insert_coord2bounds(0, nullbounds());

  gbxmaps.insert_coord3nghbrs(0, nullnghbrs(0));
  gbxmaps.insert_coord1nghbrs(0, nullnghbrs(0));
  gbxmaps.insert_coord2nghbrs(0, nullnghbrs(0));
}

/* Gives all coord[X]bounds maps for X = x or y
null values (max/min numerical limits) for all gridboxes and
neighbours maps are null (meaning periodic boundary conditions
where neighbour of gridbox in x or y direction is itself).
coord3bounds map, ie. z direction, is set using gfb.
coord3 back / forward neighbours set using finite or
periodic boundary conditions in a cartesian domain */
void set_1Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps) {
  const auto ndims(gfb.ndims);

  for (auto idx : gfb.gbxidxs) {
    const auto c3bs(gfb.get_coord3gbxbounds(idx));
    gbxmaps.insert_coord3bounds(idx, c3bs);
    gbxmaps.insert_coord1bounds(idx, nullbounds());
    gbxmaps.insert_coord2bounds(idx, nullbounds());

    const auto c3nghbrs(DoublyPeriodicDomain::cartesian_coord3nghbrs(idx, ndims));
    gbxmaps.insert_coord3nghbrs(idx, c3nghbrs);
    gbxmaps.insert_coord1nghbrs(idx, nullnghbrs(0));
    gbxmaps.insert_coord2nghbrs(idx, nullnghbrs(0));
  }
}

/* Gives coordybounds map null values (max/min numerical limits)
for all gridboxes and y neighbours maps are null (meaning periodic
boundary conditions where neighbour of gridbox in y direction is
itself). coord3 and coord1 bounds maps (ie. z and x) are set
using gfb. coord3 and coord1 neighbours maps call functions for
appropriate periodic or finite boundary conditions in a
cartesian domain  */
void set_2Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps) {
  const auto ndims(gfb.ndims);

  for (auto idx : gfb.gbxidxs) {
    const auto c3bs(gfb.get_coord3gbxbounds(idx));
    gbxmaps.insert_coord3bounds(idx, c3bs);
    const auto c1bs(gfb.get_coord1gbxbounds(idx));
    gbxmaps.insert_coord1bounds(idx, c1bs);
    gbxmaps.insert_coord2bounds(idx, nullbounds());

    const auto c3nghbrs(DoublyPeriodicDomain::cartesian_coord3nghbrs(idx, ndims));
    gbxmaps.insert_coord3nghbrs(idx, c3nghbrs);
    const auto c1nghbrs(DoublyPeriodicDomain::cartesian_coord1nghbrs(idx, ndims));
    gbxmaps.insert_coord1nghbrs(idx, c1nghbrs);
    gbxmaps.insert_coord2nghbrs(idx, nullnghbrs(0));
  }
}

kkpair_size_t correct_neighbor_indices(kkpair_size_t neighbours, const std::vector<size_t> ndims,
                                       CartesianDecomposition &domain_decomposition) {
  // If the neighbour index is not local sum total_local_gridboxes so that it
  // can be identified later If the neighbour index is local convert it to a
  // local index by subtracting gridboxes_slice_start
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  std::array<size_t, 3> neighbor_coordinates;

  if (neighbours.first != LIMITVALUES::oob_gbxindex) {
    neighbor_coordinates = get_coordinates_from_index(ndims, neighbours.first);
    if (domain_decomposition.check_indices_inside_partition(neighbor_coordinates, my_rank))
      neighbours.first = domain_decomposition.global_to_local_gridbox_index(neighbours.first);
    else
      neighbours.first += domain_decomposition.get_total_global_gridboxes();
  }

  if (neighbours.second != LIMITVALUES::oob_gbxindex) {
    neighbor_coordinates = get_coordinates_from_index(ndims, neighbours.second);
    if (domain_decomposition.check_indices_inside_partition(neighbor_coordinates, my_rank))
      neighbours.second = domain_decomposition.global_to_local_gridbox_index(neighbours.second);
    else
      neighbours.second += domain_decomposition.get_total_global_gridboxes();
  }

  return neighbours;
}

/* Sets all coord[X]bounds maps (for X = x, y, z)
using gfb data as well as back and forward neighbours
maps assuming periodic or finite boundary conditions
in cartesian domain */
void set_3Dmodel_maps(const GbxBoundsFromBinary &gfb, CartesianMaps &gbxmaps) {
  const auto ndims(gfb.ndims);
  int my_rank, local_gbx_index;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  auto domain_decomposition = gbxmaps.get_domain_decomposition();
  auto partition_origin = domain_decomposition.get_local_partition_origin();
  auto partition_size = domain_decomposition.get_local_partition_size();
  domain_decomposition.set_dimensions_bound_behavior({0, 1, 1});

  for (size_t k = 0; k < partition_size[0]; k++)
    for (size_t i = 0; i < partition_size[1]; i++)
      for (size_t j = 0; j < partition_size[2]; j++) {
        int idx = get_index_from_coordinates(ndims, partition_origin[0] + k,
                                             partition_origin[1] + i, partition_origin[2] + j);

        local_gbx_index = domain_decomposition.global_to_local_gridbox_index(idx);

        gbxmaps.insert_gbxarea(local_gbx_index, gfb.gbxarea(idx));
        gbxmaps.insert_gbxvolume(local_gbx_index, gfb.gbxvol(idx));

        const auto c3bs(gfb.get_coord3gbxbounds(idx));
        gbxmaps.insert_coord3bounds(local_gbx_index, c3bs);
        const auto c3nghbrs(DoublyPeriodicDomain::cartesian_coord3nghbrs(idx, ndims));
        gbxmaps.insert_coord3nghbrs(
            local_gbx_index, correct_neighbor_indices(c3nghbrs, ndims, domain_decomposition));

        const auto c1bs(gfb.get_coord1gbxbounds(idx));
        gbxmaps.insert_coord1bounds(local_gbx_index, c1bs);
        const auto c1nghbrs(DoublyPeriodicDomain::cartesian_coord1nghbrs(idx, ndims));
        gbxmaps.insert_coord1nghbrs(
            local_gbx_index, correct_neighbor_indices(c1nghbrs, ndims, domain_decomposition));

        const auto c2bs(gfb.get_coord2gbxbounds(idx));
        gbxmaps.insert_coord2bounds(local_gbx_index, c2bs);
        const auto c2nghbrs(DoublyPeriodicDomain::cartesian_coord2nghbrs(idx, ndims));
        gbxmaps.insert_coord2nghbrs(
            local_gbx_index, correct_neighbor_indices(c2nghbrs, ndims, domain_decomposition));
      }
}
