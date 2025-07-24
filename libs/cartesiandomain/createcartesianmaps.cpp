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

void check_ngridboxes_matches_maps(const CartesianMaps &gbxmaps, const size_t ngbxs);

void check_ngridboxes_matches_ndims(const CartesianMaps &gbxmaps, const size_t ngbxs);

void set_maps_ndims(const std::vector<size_t> &ndims, CartesianMaps &gbxmaps);

void set_cartesian_maps(const unsigned int nspacedims, const GbxBoundsFromBinary &gfb,
                        CartesianMaps &gbxmaps);

void set_null_cartesian_maps(const unsigned int nspacedims, const GbxBoundsFromBinary &gfb,
                             CartesianMaps &gbxmaps);

/* bounds for CartesianMaps of gridboxes along directions of model not used e.g. in 1-D model,
these are bounds of gridboxes in coord1 and coord2 directions */
Kokkos::pair<double, double> nullbounds() { return {LIMITVALUES::llim, LIMITVALUES::ulim}; }

/* {back, forward} neighbours for CartesianMaps of gridboxes along directions of model not used.
Boundaries are 'periodic' BCs in non-existent dimensions e.g. in 2-D model, neighbour in
coord2 direction of gridbox is itself */
Kokkos::pair<unsigned int, unsigned int> nullnghbrs(const unsigned int idx) { return {idx, idx}; }

/* creates cartesian maps instance using gridbox bounds read from gridfile for a 0-D, 1-D, 2-D
or 3-D model with periodic or finite boundary conditions. In a non-3D case, boundaries and
neighbours maps for unused dimensions are 'null' (ie. return numerical limits), however the area
and volume of each gridbox remains finite. E.g. In the 0-D case, the bounds maps all have 1
{key, value} where key=gbxidx=0 and value = {max, min} numerical limits, meanwhile volume
function returns a value determined from the gridfile 'grid_filename' */
CartesianMaps create_cartesian_maps(const size_t ngbxs, const unsigned int nspacedims,
                                    const std::filesystem::path grid_filename) {
  std::cout << "\n--- create cartesian gridbox maps ---\n";

  const auto gfb = GbxBoundsFromBinary(ngbxs, nspacedims, grid_filename);

  auto gbxmaps = CartesianMaps();

  gbxmaps.create_decomposition(gfb.ndims, gfb);
  set_cartesian_maps(nspacedims, gfb, gbxmaps);

  set_maps_ndims(gfb.ndims, gbxmaps);

  check_ngridboxes_matches_ndims(gbxmaps, gbxmaps.get_total_global_ngridboxes());
  check_ngridboxes_matches_maps(gbxmaps, gbxmaps.get_local_ngridboxes_hostcopy());

  std::cout << "--- create cartesian gridbox maps: success ---\n";

  return gbxmaps;
}

void check_ngridboxes_matches_maps(const CartesianMaps &gbxmaps, const size_t ngbxs) {
  const auto ngbxs_from_maps = gbxmaps.maps_size();
  if (ngbxs_from_maps != ngbxs + 1) {
    throw std::invalid_argument(
        "ngbxs from gridbox maps inconsistent "
        " with number of gridboxes");
  }
}

/* checks number of gridboxes according to maps matches with expected value from gfb */
void check_ngridboxes_matches_ndims(const CartesianMaps &gbxmaps, const size_t ngbxs) {
  const auto h_ndims = gbxmaps.get_global_ndims_hostcopy();
  const auto ngbxs_from_ndims = size_t{h_ndims(0) * h_ndims(1) * h_ndims(2)};

  if (ngbxs_from_ndims != ngbxs) {
    throw std::invalid_argument(
        "ndims from gridbox maps inconsistent "
        " with number of gridboxes");
  }
}

/*
copys ndims  to gbxmaps' ndims to set number of dimensions (ie. number of gridboxes) in
[coord3, coord1, coord2] directions
*/
void set_maps_ndims(const std::vector<size_t> &i_ndims, CartesianMaps &gbxmaps) {
  auto h_ndims = Kokkos::create_mirror_view(
      gbxmaps.get_global_ndims());  // mirror ndims in case view is on device

  for (unsigned int m(0); m < 3; ++m) {
    h_ndims(m) = i_ndims.at(m);
  }

  gbxmaps.set_global_ndims_via_copy(h_ndims);
}

/*
If the neighbour index is not local sum total_local_gridboxes so that it
can be identified later If the neighbour index is local convert it to a
local index by subtracting gridboxes_slice_start
*/
kkpair_size_t correct_neighbor_indices(kkpair_size_t neighbours, const std::vector<size_t> ndims,
                                       const CartesianDecomposition &domain_decomposition) {
  int my_rank;
  my_rank = init_communicator::get_comm_rank();
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

/* Sets all coord[X]bounds maps (for X = x, y, z) using gfb data as well as back and
forward neighbours maps assuming periodic or finite boundary conditions in
cartesian domain */
void set_cartesian_maps(const unsigned int nspacedims, const GbxBoundsFromBinary &gfb,
                        CartesianMaps &gbxmaps) {
  if (nspacedims > 3) {
    throw std::invalid_argument("only 0 <= nspacedims <= 3 is valid ");
  }

  const auto ndims(gfb.ndims);

  auto domain_decomposition = gbxmaps.get_domain_decomposition();
  auto partition_origin = domain_decomposition.get_local_partition_origin();
  auto partition_size = domain_decomposition.get_local_partition_size();
  domain_decomposition.set_dimensions_bound_behavior({0, 1, 1});

  const auto sz = gbxmaps.get_local_ngridboxes_hostcopy() + 1;  // +1 for oob_gbxindex key

  const auto h_to_coord3bounds = kokkos_pairmap::HostMirror(sz);
  const auto h_to_coord1bounds = kokkos_pairmap::HostMirror(sz);
  const auto h_to_coord2bounds = kokkos_pairmap::HostMirror(sz);

  const auto h_to_back_coord3nghbr = kokkos_uintmap::HostMirror(sz);
  const auto h_to_forward_coord3nghbr = kokkos_uintmap::HostMirror(sz);
  const auto h_to_back_coord1nghbr = kokkos_uintmap::HostMirror(sz);
  const auto h_to_forward_coord1nghbr = kokkos_uintmap::HostMirror(sz);
  const auto h_to_back_coord2nghbr = kokkos_uintmap::HostMirror(sz);
  const auto h_to_forward_coord2nghbr = kokkos_uintmap::HostMirror(sz);

  const auto to_gbxareas = kokkos_dblmaph(sz);
  const auto to_gbxvolumes = kokkos_dblmaph(sz);

  /* sets value for coordinate bounds, neighbours, areas and volumes
  for case when outofbounds gbxidx "oob_gbxindex" searches map */
  const auto oob = LIMITVALUES::oob_gbxindex;
  h_to_coord3bounds.insert(oob, nullbounds());
  h_to_coord1bounds.insert(oob, nullbounds());
  h_to_coord2bounds.insert(oob, nullbounds());
  h_to_back_coord3nghbr.insert(oob, nullnghbrs(oob).first);
  h_to_forward_coord3nghbr.insert(oob, nullnghbrs(oob).second);
  h_to_back_coord1nghbr.insert(oob, nullnghbrs(oob).first);
  h_to_forward_coord1nghbr.insert(oob, nullnghbrs(oob).second);
  h_to_back_coord2nghbr.insert(oob, nullnghbrs(oob).first);
  h_to_forward_coord2nghbr.insert(oob, nullnghbrs(oob).second);
  to_gbxareas.insert(oob, 0.0);
  to_gbxvolumes.insert(oob, 0.0);

  Kokkos::parallel_for(
      "set_cartesian_maps", HostTeamPolicy(partition_size[0], Kokkos::AUTO()),
      KOKKOS_LAMBDA(const HostTeamMember &team_member) {
        const auto k = team_member.league_rank();
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_member, partition_size[1]), [=](const size_t i) {
              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(team_member, partition_size[2]), [=](const size_t j) {
                    int idx = get_index_from_coordinates(ndims, partition_origin[0] + k,
                                                         partition_origin[1] + i,
                                                         partition_origin[2] + j);

                    int local_gbx_index = domain_decomposition.global_to_local_gridbox_index(idx);

                    const auto c3bs = gfb.get_coord3gbxbounds(idx);
                    h_to_coord3bounds.insert(local_gbx_index, c3bs);
                    auto c3nghbrs = DoublyPeriodicDomain::cartesian_coord3nghbrs(idx, ndims);
                    c3nghbrs = correct_neighbor_indices(c3nghbrs, ndims, domain_decomposition);
                    h_to_back_coord3nghbr.insert(local_gbx_index, c3nghbrs.first);
                    h_to_forward_coord3nghbr.insert(local_gbx_index, c3nghbrs.second);

                    const auto c1bs = gfb.get_coord1gbxbounds(idx);
                    h_to_coord1bounds.insert(local_gbx_index, c1bs);
                    auto c1nghbrs = DoublyPeriodicDomain::cartesian_coord1nghbrs(idx, ndims);
                    c1nghbrs = correct_neighbor_indices(c1nghbrs, ndims, domain_decomposition);
                    h_to_back_coord1nghbr.insert(local_gbx_index, c1nghbrs.first);
                    h_to_forward_coord1nghbr.insert(local_gbx_index, c1nghbrs.second);

                    const auto c2bs = gfb.get_coord2gbxbounds(idx);
                    h_to_coord2bounds.insert(local_gbx_index, c2bs);
                    auto c2nghbrs = DoublyPeriodicDomain::cartesian_coord2nghbrs(idx, ndims);
                    c2nghbrs = correct_neighbor_indices(c2nghbrs, ndims, domain_decomposition);
                    h_to_back_coord2nghbr.insert(local_gbx_index, c2nghbrs.first);
                    h_to_forward_coord2nghbr.insert(local_gbx_index, c2nghbrs.second);

                    to_gbxareas.insert(local_gbx_index, gfb.gbxarea(idx));
                    to_gbxvolumes.insert(local_gbx_index, gfb.gbxvol(idx));
                  });
            });
      });

  switch (nspacedims) {
    case 3:  // 3-D model (set coord2 dimension)
      gbxmaps.set_coord2bounds_via_copy(h_to_coord2bounds);
      gbxmaps.set_back_coord2nghbr_via_copy(h_to_back_coord2nghbr);
      gbxmaps.set_forward_coord2nghbr_via_copy(h_to_forward_coord2nghbr);
      [[fallthrough]];
    case 2:  // 3-D or 2-D model (set coord1 dimension)
      gbxmaps.set_coord1bounds_via_copy(h_to_coord1bounds);
      gbxmaps.set_back_coord1nghbr_via_copy(h_to_back_coord1nghbr);
      gbxmaps.set_forward_coord1nghbr_via_copy(h_to_forward_coord1nghbr);
      [[fallthrough]];
    case 1:  // 3-D, 2-D or 1-D model (set coord3 dimension)
      gbxmaps.set_coord3bounds_via_copy(h_to_coord3bounds);
      gbxmaps.set_back_coord3nghbr_via_copy(h_to_back_coord3nghbr);
      gbxmaps.set_forward_coord3nghbr_via_copy(h_to_forward_coord3nghbr);
      [[fallthrough]];
    case 0:  // 3-D, 2-D, 1-D or 0-D model (set areas and volumes)
      gbxmaps.set_gbxareas_map(to_gbxareas);
      gbxmaps.set_gbxvolumes_map(to_gbxvolumes);
  }

  if (nspacedims < 3) {
    set_null_cartesian_maps(nspacedims, gfb, gbxmaps);
  }
}

/*
For null dimensions (see below), function gives coord[X]bounds maps null values (max/min numerical
limits) for all gridboxes and also gives neighbours maps null values (meaning periodic boundary
conditions where neighbour of gridbox in a certain direction is itself). Null dimensions are:
 - coord2 (y) for a 2-D model,
 - coord1 and coord2 (x and y) for a 1-D model,
 - coord3, coord1 and coord2 (z, x and y) for a 0-D model.
*/
void set_null_cartesian_maps(const unsigned int nspacedims, const GbxBoundsFromBinary &gfb,
                             CartesianMaps &gbxmaps) {
  if (nspacedims >= 3) {
    throw std::invalid_argument("null model dimensions only valid for 0 <= nspacedims < 3");
  }

  const auto ndims(gfb.ndims);

  auto domain_decomposition = gbxmaps.get_domain_decomposition();
  auto partition_origin = domain_decomposition.get_local_partition_origin();
  auto partition_size = domain_decomposition.get_local_partition_size();
  domain_decomposition.set_dimensions_bound_behavior({0, 1, 1});

  const auto sz = gbxmaps.get_local_ngridboxes_hostcopy() + 1;  // +1 for oob_gbxindex key
  const auto h_nullbounds = kokkos_pairmap::HostMirror(sz);
  const auto h_back_nullnghbr = kokkos_uintmap::HostMirror(sz);
  const auto h_forward_nullnghbr = kokkos_uintmap::HostMirror(sz);

  /* sets value for coordinate bounds, neighbours, areas and volumes
  for case when outofbounds gbxidx "oob_gbxindex" searches map */
  const auto oob = LIMITVALUES::oob_gbxindex;
  h_nullbounds.insert(oob, nullbounds());
  h_back_nullnghbr.insert(oob, nullnghbrs(oob).first);
  h_forward_nullnghbr.insert(oob, nullnghbrs(oob).first);

  // TODO(ALL): perform on GPUs once domain_decomposition is GPU compatible (then after loop assign
  // gbxmaps in switch rather than use deep_copy)
  Kokkos::parallel_for(
      "set_null_cartesian_maps", HostTeamPolicy(partition_size[0], Kokkos::AUTO()),
      KOKKOS_LAMBDA(const HostTeamMember &team_member) {
        const auto k = team_member.league_rank();
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_member, partition_size[1]), [=](const size_t i) {
              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(team_member, partition_size[2]), [=](const size_t j) {
                    int idx = get_index_from_coordinates(ndims, partition_origin[0] + k,
                                                         partition_origin[1] + i,
                                                         partition_origin[2] + j);

                    int local_gbx_index = domain_decomposition.global_to_local_gridbox_index(idx);

                    h_nullbounds.insert(local_gbx_index, nullbounds());
                    const auto nullnghbr = nullnghbrs(local_gbx_index);
                    h_back_nullnghbr.insert(local_gbx_index, nullnghbr.first);
                    h_forward_nullnghbr.insert(local_gbx_index, nullnghbr.second);
                  });
            });
      });

  switch (nspacedims) {
    case 0:  // 0-D model (set coord3 dimension null)
      gbxmaps.set_coord3bounds_via_copy(h_nullbounds);
      gbxmaps.set_back_coord3nghbr_via_copy(h_back_nullnghbr);
      gbxmaps.set_forward_coord3nghbr_via_copy(h_forward_nullnghbr);
      [[fallthrough]];
    case 1:  // 1-D or 0-D model (set coord1 dimension null)
      gbxmaps.set_coord1bounds_via_copy(h_nullbounds);
      gbxmaps.set_back_coord1nghbr_via_copy(h_back_nullnghbr);
      gbxmaps.set_forward_coord1nghbr_via_copy(h_forward_nullnghbr);
      [[fallthrough]];
    case 2:  // 2-D, 1-D or 0-D model (set coord2 dimension null)
      gbxmaps.set_coord2bounds_via_copy(h_nullbounds);
      gbxmaps.set_back_coord2nghbr_via_copy(h_back_nullnghbr);
      gbxmaps.set_forward_coord2nghbr_via_copy(h_forward_nullnghbr);
  }
}
