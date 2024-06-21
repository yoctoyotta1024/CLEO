/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesianmotion.cpp
 * Project: cartesiandomain
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 */

#include "./cartesianmotion.hpp"

KOKKOS_FUNCTION int flag_sdgbxindex(const unsigned int idx,
                                    const Kokkos::pair<double, double> bounds, const double coord);

KOKKOS_FUNCTION unsigned int change_to_backwards_coord3nghbr(const unsigned int idx,
                                                             const CartesianMaps &gbxmaps,
                                                             Superdrop &superdrop);
KOKKOS_FUNCTION void beyonddomain_backwards_coord3(const CartesianMaps &gbxmaps,
                                                   const unsigned int idx, const unsigned int nghbr,
                                                   Superdrop &drop);

KOKKOS_FUNCTION unsigned int change_to_forwards_coord3nghbr(const unsigned int idx,
                                                            const CartesianMaps &gbxmaps,
                                                            Superdrop &drop);
KOKKOS_FUNCTION void beyonddomain_forwards_coord3(const CartesianMaps &gbxmaps,
                                                  const unsigned int idx, const unsigned int nghbr,
                                                  Superdrop &drop);

KOKKOS_FUNCTION unsigned int change_to_backwards_coord1nghbr(const unsigned int idx,
                                                             const CartesianMaps &gbxmaps,
                                                             Superdrop &superdrop);
KOKKOS_FUNCTION void beyonddomain_backwards_coord1(const CartesianMaps &gbxmaps,
                                                   const unsigned int idx, const unsigned int nghbr,
                                                   Superdrop &drop);

KOKKOS_FUNCTION unsigned int change_to_forwards_coord1nghbr(const unsigned int idx,
                                                            const CartesianMaps &gbxmaps,
                                                            Superdrop &drop);
KOKKOS_FUNCTION void beyonddomain_forwards_coord1(const CartesianMaps &gbxmaps,
                                                  const unsigned int idx, const unsigned int nghbr,
                                                  Superdrop &drop);

KOKKOS_FUNCTION unsigned int change_to_backwards_coord2nghbr(const unsigned int idx,
                                                             const CartesianMaps &gbxmaps,
                                                             Superdrop &superdrop);
KOKKOS_FUNCTION void beyonddomain_backwards_coord2(const CartesianMaps &gbxmaps,
                                                   const unsigned int idx, const unsigned int nghbr,
                                                   Superdrop &drop);

KOKKOS_FUNCTION unsigned int change_to_forwards_coord2nghbr(const unsigned int idx,
                                                            const CartesianMaps &gbxmaps,
                                                            Superdrop &drop);
KOKKOS_FUNCTION void beyonddomain_forwards_coord2(const CartesianMaps &gbxmaps,
                                                  const unsigned int idx, const unsigned int nghbr,
                                                  Superdrop &drop);

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
int flag_sdgbxindex(const unsigned int idx, const Kokkos::pair<double, double> bounds,
                    const double coord) {
  if (idx == outofbounds_gbxindex()) {
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

/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord3 direction.
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
_Note:_ backwards/forwards functions may change the
superdroplet's attributes e.g. if it leaves the domain. */
KOKKOS_FUNCTION unsigned int change_if_coord3nghbr(const CartesianMaps &gbxmaps, unsigned int idx,
                                                   Superdrop &drop) {
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord3bounds(idx),
                                    drop.get_coord3());  // if value != 0 idx needs to change
  switch (flag) {
    case 1:
      idx = change_to_backwards_coord3nghbr(idx, gbxmaps, drop);
      break;
    case 2:
      idx = change_to_forwards_coord3nghbr(idx, gbxmaps, drop);
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
KOKKOS_FUNCTION unsigned int change_if_coord1nghbr(const CartesianMaps &gbxmaps, unsigned int idx,
                                                   Superdrop &drop) {
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord1bounds(idx),
                                    drop.get_coord1());  // if value != 0 idx needs to change
  switch (flag) {
    case 1:
      idx = change_to_backwards_coord1nghbr(idx, gbxmaps, drop);
      break;
    case 2:
      idx = change_to_forwards_coord1nghbr(idx, gbxmaps, drop);
      break;
  }
  return idx;
}

/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord2 direction.
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
_Note:_ backwards/forwards functions may change the
superdroplet's attributes e.g. if it leaves the domain. */
KOKKOS_FUNCTION unsigned int change_if_coord2nghbr(const CartesianMaps &gbxmaps, unsigned int idx,
                                                   Superdrop &drop) {
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord2bounds(idx),
                                    drop.get_coord2());  // if value != 0 idx needs to change
  switch (flag) {
    case 1:
      idx = change_to_backwards_coord2nghbr(idx, gbxmaps, drop);
      break;
    case 2:
      idx = change_to_forwards_coord2nghbr(idx, gbxmaps, drop);
      break;
  }
  return idx;
}

/* function to return gbxindex of neighbouring gridbox
in backwards coord3 (z) direction and to update superdrop
if its coord3 has exceeded the z lower domain boundary */
KOKKOS_FUNCTION unsigned int change_to_backwards_coord3nghbr(const unsigned int idx,
                                                             const CartesianMaps &gbxmaps,
                                                             Superdrop &drop) {
  const auto nghbr = (unsigned int)gbxmaps.coord3backward(idx);

  const auto incre = (unsigned int)1;  // increment
  // drop was at lower z edge of domain (now moving below it)
  if (beyond_domainboundary(idx, incre, gbxmaps.get_ndim(0))) {
    beyonddomain_backwards_coord3(gbxmaps, idx, nghbr, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of z backwards (down) neighbour
}

/* function to return gbxindex of neighbouring gridbox in
forwards coord3 (z) direction and to update superdrop coord3
if superdrop has exceeded the z upper domain boundary */
KOKKOS_FUNCTION unsigned int change_to_forwards_coord3nghbr(const unsigned int idx,
                                                            const CartesianMaps &gbxmaps,
                                                            Superdrop &drop) {
  const auto nghbr = (unsigned int)gbxmaps.coord3forward(idx);

  const auto incre = (unsigned int)1;  // increment
  // drop was upper z edge of domain (now moving above it)
  if (beyond_domainboundary(idx + incre, incre, gbxmaps.get_ndim(0))) {
    beyonddomain_forwards_coord3(gbxmaps, idx, nghbr, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of z forwards (up) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in backwards coord1 (x) direction and to update superdrop
coord1 if superdrop has exceeded the x back domain boundary */
KOKKOS_FUNCTION unsigned int change_to_backwards_coord1nghbr(const unsigned int idx,
                                                             const CartesianMaps &gbxmaps,
                                                             Superdrop &drop) {
  const auto nghbr = (unsigned int)gbxmaps.coord1backward(idx);

  const auto ndims(gbxmaps.get_ndims());
  const auto incre = (unsigned int)ndims(0);  // increment
  // at lower x edge of domain
  if (beyond_domainboundary(idx, incre, ndims(1))) {
    beyonddomain_backwards_coord1(gbxmaps, idx, nghbr, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of x backwards (behind) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in forwards coord1 (x) direction and to update superdrop
coord1 if superdrop has exceeded the x front domain boundary */
KOKKOS_FUNCTION unsigned int change_to_forwards_coord1nghbr(const unsigned int idx,
                                                            const CartesianMaps &gbxmaps,
                                                            Superdrop &drop) {
  const auto nghbr = (unsigned int)gbxmaps.coord1forward(idx);

  const auto ndims(gbxmaps.get_ndims());
  const auto incre = (unsigned int)ndims(0);  // increment
  // at lower x edge of domain
  if (beyond_domainboundary(idx + incre, incre, ndims(1))) {
    beyonddomain_forwards_coord1(gbxmaps, idx, nghbr, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of x forwards (infront) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in backwards coord2 (y) direction and to update superdrop
coord2 if superdrop has exceeded the y leftmost domain boundary */
KOKKOS_FUNCTION unsigned int change_to_backwards_coord2nghbr(const unsigned int idx,
                                                             const CartesianMaps &gbxmaps,
                                                             Superdrop &drop) {
  const auto nghbr = (unsigned int)gbxmaps.coord2backward(idx);

  const auto ndims(gbxmaps.get_ndims());
  const auto incre = (unsigned int)ndims(0) *
                     ndims(1);  // no. gridboxes in z direction * no. gridboxes in x direction
  // at lower y edge of domain
  if (beyond_domainboundary(idx, incre, ndims(2))) {
    beyonddomain_backwards_coord2(gbxmaps, idx, nghbr, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of y backwards (left) neighbour
}

/* function to return gbxindex of neighbouring gridbox
in forwards coord2 (y) direction and to update superdrop
coord2 if superdrop has exceeded the y rightmost domain boundary */
KOKKOS_FUNCTION unsigned int change_to_forwards_coord2nghbr(const unsigned int idx,
                                                            const CartesianMaps &gbxmaps,
                                                            Superdrop &drop) {
  const auto nghbr = (unsigned int)gbxmaps.coord2forward(idx);

  const auto ndims(gbxmaps.get_ndims());
  const auto incre = (unsigned int)ndims(0) *
                     ndims(1);  // no. gridboxes in z direction * no. gridboxes in x direction
  // at upper y edge of domain
  if (beyond_domainboundary(idx + incre, incre, ndims(2))) {
    beyonddomain_forwards_coord2(gbxmaps, idx, nghbr, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of y forwards (right) neighbour
}

/* function updates superdrop that has crossed the
domain boundary in the backwards coord3 (z) direction
(i.e. superdrop has exceeded the z lower domain boundary) */
KOKKOS_FUNCTION void beyonddomain_backwards_coord3(const CartesianMaps &gbxmaps,
                                                   const unsigned int idx, const unsigned int nghbr,
                                                   Superdrop &drop) {
  const auto lim1 = double{gbxmaps.coord3bounds(nghbr).second};  // upper lim of backward neighbour
  const auto lim2 = double{gbxmaps.coord3bounds(idx).first};     // lower lim of current gbx
  drop.set_coord3(DoublyPeriodicDomain::boundarycond_coord3(drop.get_coord3(), lim1, lim2));
}

/* function updates superdrop that has crossed the
domain boundary in the forwards coord3 (z) direction
(i.e. superdrop has exceeded the z upper domain boundary) */
KOKKOS_FUNCTION void beyonddomain_forwards_coord3(const CartesianMaps &gbxmaps,
                                                  const unsigned int idx, const unsigned int nghbr,
                                                  Superdrop &drop) {
  const auto lim1 = double{gbxmaps.coord3bounds(nghbr).first};  // lower lim of forward neighbour
  const auto lim2 = double{gbxmaps.coord3bounds(idx).second};   // upper lim of current gbx
  drop.set_coord3(DoublyPeriodicDomain::boundarycond_coord3(drop.get_coord3(), lim1, lim2));
}

/* function updates superdrop that has crossed the
domain boundary in the backwards coord1 (x) direction
(i.e. superdrop has exceeded the x back domain boundary) */
KOKKOS_FUNCTION void beyonddomain_backwards_coord1(const CartesianMaps &gbxmaps,
                                                   const unsigned int idx, const unsigned int nghbr,
                                                   Superdrop &drop) {
  const auto lim1 = double{gbxmaps.coord1bounds(nghbr).second};  // upper lim of backward neigghbour
  const auto lim2 = double{gbxmaps.coord1bounds(idx).first};     // lower lim of current gbx
  drop.set_coord1(DoublyPeriodicDomain::boundarycond_coord1(drop.get_coord1(), lim1, lim2));
}

/* function updates superdrop that has crossed the
domain boundary in the forwards coord1 (x) direction
(i.e. superdrop has exceeded the x front domain boundary) */
KOKKOS_FUNCTION void beyonddomain_forwards_coord1(const CartesianMaps &gbxmaps,
                                                  const unsigned int idx, const unsigned int nghbr,
                                                  Superdrop &drop) {
  const auto lim1 = double{gbxmaps.coord1bounds(nghbr).first};  // lower lim of forward nghbour
  const auto lim2 = double{gbxmaps.coord1bounds(idx).second};   // upper lim of gbx
  drop.set_coord1(DoublyPeriodicDomain::boundarycond_coord1(drop.get_coord1(), lim1, lim2));
}

/* function updates superdrop that has crossed the
domain boundary in the backwards coord2 (y) direction
(i.e. superdrop has exceeded the y leftmost domain boundary) */
KOKKOS_FUNCTION void beyonddomain_backwards_coord2(const CartesianMaps &gbxmaps,
                                                   const unsigned int idx, const unsigned int nghbr,
                                                   Superdrop &drop) {
  const auto lim1 = double{gbxmaps.coord2bounds(nghbr).second};  // upper lim of backward nghbour
  const auto lim2 = double{gbxmaps.coord2bounds(idx).first};     // lower lim of gbx
  drop.set_coord2(DoublyPeriodicDomain::boundarycond_coord2(drop.get_coord2(), lim1, lim2));
}

/* function updates superdrop that has crossed the
domain boundary in the forwards coord2 (y) direction
(i.e. superdrop has exceeded the y rightmost domain boundary) */
KOKKOS_FUNCTION void beyonddomain_forwards_coord2(const CartesianMaps &gbxmaps,
                                                  const unsigned int idx, const unsigned int nghbr,
                                                  Superdrop &drop) {
  const auto lim1 = double{gbxmaps.coord2bounds(nghbr).first};  // lower lim of forward nghbour
  const auto lim2 = double{gbxmaps.coord2bounds(idx).second};   // upper lim of gbx
  drop.set_coord2(DoublyPeriodicDomain::boundarycond_coord2(drop.get_coord2(), lim1, lim2));
}
