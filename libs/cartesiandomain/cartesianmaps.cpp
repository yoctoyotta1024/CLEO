/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesianmaps.cpp
 * Project: cartesiandomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions related to creating and using maps to convert
 * between a gridbox indexes and domain coordinates for a
 * cartesian C grid
 */

#include "cartesiandomain/cartesianmaps.hpp"

KOKKOS_FUNCTION
unsigned int get_no_decomposition_bounding_gridbox(const CartesianMaps &gbxmaps,
                                                   const unsigned int gbxindex, double &coord3,
                                                   double &coord1, double &coord2);

/* on host, throws error if maps are not all
the same size, else returns size of maps */
size_t CartesianMaps::maps_size() const {
  // ngbxs + 1 for out of bounds key
  const auto sz = domain_decomposition.get_total_local_gridboxes() + 1;

  if (to_coord3bounds.size() != sz || to_coord1bounds.size() != sz ||
      to_coord2bounds.size() != sz || to_back_coord3nghbr.size() != sz ||
      to_forward_coord3nghbr.size() != sz || to_back_coord1nghbr.size() != sz ||
      to_forward_coord1nghbr.size() != sz || to_back_coord2nghbr.size() != sz ||
      to_forward_coord2nghbr.size() != sz || to_areas.size() != sz || to_volumes.size() != sz) {
    throw std::invalid_argument("gridbox maps are not all the same size");
  }

  return sz;
}

// TODO(ALL) make domain_decomp call compatible with GPUs and then remove comm_size guard
KOKKOS_FUNCTION
size_t CartesianMaps::get_local_ngridboxes() const {
  if (is_decomp) {
    return domain_decomposition.get_total_local_gridboxes();
  }
  return global_ndims(0) * global_ndims(1) * global_ndims(2);
}

// TODO(ALL) make domain_decomp call compatible with GPUs and then remove comm_size guard
KOKKOS_FUNCTION
size_t CartesianMaps::local_to_global_gridbox_index(unsigned int local_gridbox_index,
                                                    int process) const {
  if (is_decomp) {
    return domain_decomposition.local_to_global_gridbox_index(local_gridbox_index, process);
  }
  return local_gridbox_index;
}

/* given coordinates, associated gxbindex is found. The coords may be updated too,
 * e.g. if the domain has a cyclic boundary condition and they therefore need to be corrected
 */
// TODO(ALL) make domain_decomp call compatible with GPUs and then remove comm_size guard
KOKKOS_FUNCTION
unsigned int CartesianMaps::get_local_bounding_gridbox_index(const unsigned int gbxindex,
                            double &coord3, double &coord1, double &coord2) const {
  if (is_decomp) {
    auto coordinates = std::array<double, 3>{coord3, coord1, coord2};
    const auto idx = domain_decomposition.get_local_bounding_gridbox_index(coordinates);
    coord3 = coordinates[0];
    coord1 = coordinates[1];
    coord2 = coordinates[2];
    return idx;
  }
  return get_no_decomposition_bounding_gridbox(*this, gbxindex, coord3, coord1, coord2);
}

/* returns flag to keep idx the same (flag = 0) or
update to forwards (flag = 1) or backwards (flag = 2)
neighbour. Flag = 0 if idx is out of domain value or
if coord lies within bounds = {lowerbound, upperbound}.
(_Note:_ lower bound inclusive and upper bound exclusive,
ie. lowerbound <= coord < upperbound).
Flag = 1 if coord < lowerbound, indicating idx should
be updated to backwards neighbour.
Flag = 2 if coord >= upperbound, indicating idx should
be updated to forwards neighbour. */
KOKKOS_FUNCTION
int flag_sdgbxindex(const unsigned int idx, const Kokkos::pair<double, double> bounds,
                    const double coord) {
  if (idx == LIMITVALUES::oob_gbxindex) {
    return 0;  // maintian idx that is already out of domain
    // lowerbound
  } else if (coord < bounds.first) {
    return 1;  // idx -> backwards_neighbour
    // upperbound
  } else if (coord >= bounds.second) {
    return 2;  // idx -> forwards_neighbour
  } else {
    return 0;  // maintain idx if coord within bounds
  }
}

/* function to return gbxindex of neighbouring gridbox
in forwards coord2 (y) direction and to update superdrop
coord2 if superdrop has exceeded the y rightmost domain boundary */
KOKKOS_FUNCTION
unsigned int change_to_forwards_coord2nghbr(const unsigned int idx, const CartesianMaps &gbxmaps,
                                            double &coord2) {
  const auto nghbr = gbxmaps.coord2forward(idx);
  const auto ndims = gbxmaps.get_global_ndims();
  const auto incre = (unsigned int)ndims(0) * ndims(1);  // ngbxs in z * ngbxs in x direction
  // at upper y edge of domain
  if (beyond_domainboundary(idx + incre, incre, ndims(2))) {
    const auto lim1 = gbxmaps.coord2bounds(nghbr).first;  // lower lim of forward nghbour
    const auto lim2 = gbxmaps.coord2bounds(idx).second;   // upper lim of gbx
    coord2 = DoublyPeriodicDomain::boundarycond_coord2(coord2, lim1, lim2);
  }
  return nghbr;  // gbxindex of y forwards (right) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in backwards coord2 (y) direction and to update superdrop
coord2 if superdrop has exceeded the y leftmost domain boundary */
KOKKOS_FUNCTION
unsigned int change_to_backwards_coord2nghbr(const unsigned int idx, const CartesianMaps &gbxmaps,
                                             double &coord2) {
  const auto nghbr = gbxmaps.coord2backward(idx);
  const auto ndims = gbxmaps.get_global_ndims();
  const auto incre = (unsigned int)ndims(0) * ndims(1);  // ngbxs in z * ngbxs in x direction
  // at lower y edge of domain
  if (beyond_domainboundary(idx, incre, ndims(2))) {
    const auto lim1 = gbxmaps.coord2bounds(nghbr).second;  // upper lim of backward nghbour
    const auto lim2 = gbxmaps.coord2bounds(idx).first;     // lower lim of gbx
    coord2 = DoublyPeriodicDomain::boundarycond_coord2(coord2, lim1, lim2);
  }
  return nghbr;  // gbxindex of y backwards (left) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in forwards coord1 (x) direction and to update superdrop
coord1 if superdrop has exceeded the x front domain boundary */
KOKKOS_FUNCTION
unsigned int change_to_forwards_coord1nghbr(const unsigned int idx, const CartesianMaps &gbxmaps,
                                            double &coord1) {
  const auto nghbr = gbxmaps.coord1forward(idx);
  const auto ndims = gbxmaps.get_global_ndims();
  const auto incre = (unsigned int)ndims(0);  // increment
  // at lower x edge of domain
  if (beyond_domainboundary(idx + incre, incre, ndims(1))) {
    const auto lim1 = gbxmaps.coord1bounds(nghbr).first;  // lower lim of forward nghbour
    const auto lim2 = gbxmaps.coord1bounds(idx).second;   // upper lim of gbx
    coord1 = DoublyPeriodicDomain::boundarycond_coord1(coord1, lim1, lim2);
  }
  return nghbr;  // gbxindex of x forwards (infront) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in backwards coord1 (x) direction and to update superdrop
coord1 if superdrop has exceeded the x back domain boundary */
KOKKOS_FUNCTION
unsigned int change_to_backwards_coord1nghbr(const unsigned int idx, const CartesianMaps &gbxmaps,
                                             double &coord1) {
  const auto nghbr = gbxmaps.coord1backward(idx);
  const auto ndims = gbxmaps.get_global_ndims();
  const auto incre = (unsigned int)ndims(0);  // increment
  // at lower x edge of domain
  if (beyond_domainboundary(idx, incre, ndims(1))) {
    const auto lim1 = gbxmaps.coord1bounds(nghbr).second;  // upper lim of backward neigghbour
    const auto lim2 = gbxmaps.coord1bounds(idx).first;     // lower lim of current gbx
    coord1 = DoublyPeriodicDomain::boundarycond_coord1(coord1, lim1, lim2);
  }
  return nghbr;  // gbxindex of x backwards (behind) neighbour
}

/* function to return gbxindex of neighbouring gridbox in
forwards coord3 (z) direction and to update superdrop coord3
if it has exceeded the z upper domain boundary */
KOKKOS_FUNCTION
unsigned int change_to_forwards_coord3nghbr(const unsigned int idx, const CartesianMaps &gbxmaps,
                                            double &coord3) {
  const auto nghbr = (unsigned int)gbxmaps.coord3forward(idx);
  const auto incre = (unsigned int)1;  // increment
  // drop was upper z edge of domain (now moving above it)
  if (beyond_domainboundary(idx + incre, incre, gbxmaps.get_global_ndim(0))) {
    const auto lim1 = gbxmaps.coord3bounds(nghbr).first;  // lower lim of forward neighbour
    const auto lim2 = gbxmaps.coord3bounds(idx).second;   // upper lim of current gbx
    coord3 = DoublyPeriodicDomain::boundarycond_coord3(coord3, lim1, lim2);
  }
  return nghbr;  // gbxindex of z forwards (up) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in backwards coord3 (z) direction and to update superdrop
coord3 if its has exceeded the z lower domain boundary */
KOKKOS_FUNCTION
unsigned int change_to_backwards_coord3nghbr(const unsigned int idx, const CartesianMaps &gbxmaps,
                                             double &coord3) {
  const auto nghbr = gbxmaps.coord3backward(idx);
  const auto incre = (unsigned int)1;  // increment
  // drop was at lower z edge of domain (now moving below it)
  if (beyond_domainboundary(idx, incre, gbxmaps.get_global_ndim(0))) {
    const auto lim1 = gbxmaps.coord3bounds(nghbr).second;  // upper lim of backward neighbour
    const auto lim2 = gbxmaps.coord3bounds(idx).first;     // lower lim of current gbx
    coord3 = DoublyPeriodicDomain::boundarycond_coord3(coord3, lim1, lim2);
  }
  return nghbr;  // gbxindex of z backwards (down) neighbour
}

/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord2 direction.
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
_Note:_ backwards/forwards functions may change the
superdroplet's attributes e.g. if it leaves the domain. */
KOKKOS_FUNCTION
unsigned int change_if_coord2nghbr(const CartesianMaps &gbxmaps, unsigned int idx, double &coord2) {
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord2bounds(idx), coord2);  // !=0 idx shld change
  switch (flag) {
    case 1:
      idx = change_to_backwards_coord2nghbr(idx, gbxmaps, coord2);
      break;
    case 2:
      idx = change_to_forwards_coord2nghbr(idx, gbxmaps, coord2);
      break;
  }
  return idx;
}

/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord1 direction.
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
_Note:_ backwards/forwards functions may change the
superdroplet's attributes e.g. if it leaves the domain. */
KOKKOS_FUNCTION
unsigned int change_if_coord1nghbr(const CartesianMaps &gbxmaps, unsigned int idx, double &coord1) {
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord1bounds(idx), coord1);  // !=0 idx shld change
  switch (flag) {
    case 1:
      idx = change_to_backwards_coord1nghbr(idx, gbxmaps, coord1);
      break;
    case 2:
      idx = change_to_forwards_coord1nghbr(idx, gbxmaps, coord1);
      break;
  }
  return idx;
}

/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord3 direction.
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
_Note:_ backwards/forwards functions may change the
superdroplet's coords e.g. if it leaves the domain. */
KOKKOS_FUNCTION
unsigned int change_if_coord3nghbr(const CartesianMaps &gbxmaps, unsigned int idx, double &coord3) {
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord3bounds(idx), coord3);  // !=0 idx shld change
  switch (flag) {
    case 1:
      idx = change_to_backwards_coord3nghbr(idx, gbxmaps, coord3);
      break;
    case 2:
      idx = change_to_forwards_coord3nghbr(idx, gbxmaps, coord3);
      break;
  }
  return idx;
}

KOKKOS_FUNCTION
unsigned int get_no_decomposition_bounding_gridbox(const CartesianMaps &gbxmaps,
                                                   const unsigned int gbxindex, double &coord3,
                                                   double &coord1, double &coord2) {
  auto idx = (unsigned int)change_if_coord3nghbr(gbxmaps, gbxindex, coord3);
  idx = change_if_coord1nghbr(gbxmaps, idx, coord1);
  idx = change_if_coord2nghbr(gbxmaps, idx, coord2);
  return idx;
}
