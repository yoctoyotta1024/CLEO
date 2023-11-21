/*
 * ----- CLEO -----
 * File: cartesianmotion.hpp
 * Project: cartesiandomain
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 21st November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * the sdgbxindex updated accordingly for a
 * cartesian domain with finite/periodi boundary
 * conditions
 */

#ifndef CARTESIANMOTION_HPP
#define CARTESIANMOTION_HPP

#include <functional>
#include <cassert>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "./cartesianmaps.hpp"
#include "./cartesianboundaryconds.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"
#include "gridboxes/predcorr.hpp"

KOKKOS_FUNCTION unsigned int
update_if_coord3nghbr(const CartesianMaps &gbxmaps,
                      unsigned int idx,
                      Superdrop &drop);
KOKKOS_FUNCTION unsigned int
update_if_coord1nghbr(const CartesianMaps &gbxmaps,
                      unsigned int idx,
                      Superdrop &drop);
KOKKOS_FUNCTION unsigned int
update_if_coord2nghbr(const CartesianMaps &gbxmaps,
                      unsigned int idx,
                      Superdrop &drop);

KOKKOS_FUNCTION void
check_inbounds_or_outdomain(const unsigned int idx,
                            const Kokkos::pair<double, double> bounds,
                            const double coord);

template <VelocityFormula TV>
struct CartesianMotion
/* satisfies motion concept for motion of a superdroplet
using a predictor-corrector method to update a superdroplet's
coordinates and then updating it's sdgbxindex using the
UpdateSdgbxindex struct for a cartesian domain */
{
  const unsigned int interval; // integer timestep for movement
  PredCorrMotion<CartesianMaps, TV> update_superdrop_coords;

  CartesianMotion(const unsigned int motionstep,
                  const std::function<double(unsigned int)> int2time,
                  const TV i_terminalv)
      : interval(motionstep),
        update_superdrop_coords(interval, int2time, i_terminalv) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const
  {
    return ((t_sdm / interval) + 1) * interval;
  }

  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int t_sdm) const
  {
    return t_sdm % interval == 0;
  }

  KOKKOS_INLINE_FUNCTION void
  update_superdrop_gbxindex(const unsigned int gbxindex,
                            const CartesianMaps &gbxmaps,
                            Superdrop &drop) const
  /* function satisfies requirements of
  "update_superdrop_gbxindex" in the motion concept to update a
  superdroplet if it should move between gridboxes in a
  cartesian domain. For each direction (z, then x, then y),
  superdrop coord is compared to gridbox bounds given by gbxmaps
  for the current gbxindex 'idx'. If superdrop coord lies outside
  bounds, forward or backward neighbour functions are called as
  appropriate  to update sdgbxindex (and possibly other superdrop
  attributes) */
  {
    unsigned int idx(gbxindex);

    idx = update_if_coord3nghbr(gbxmaps, idx, drop);
    check_inbounds_or_outdomain(idx, gbxmaps.coord3bounds(idx),
                                drop.get_coord3());

    idx = update_if_coord1nghbr(gbxmaps, idx, drop);
    check_inbounds_or_outdomain(idx, gbxmaps.coord1bounds(idx),
                                drop.get_coord1());

    idx = update_if_coord2nghbr(gbxmaps, idx, drop);
    check_inbounds_or_outdomain(idx, gbxmaps.coord2bounds(idx),
                                drop.get_coord2());

    drop.set_sdgbxindex(idx);
  }
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION int
flag_sdgbxindex(const unsigned int idx,
                const Kokkos::pair<double, double> bounds,
                const double coord);

KOKKOS_FUNCTION unsigned int
backwards_coord3idx(const unsigned int idx,
                    const CartesianMaps &gbxmaps,
                    Superdrop &superdrop);
KOKKOS_FUNCTION unsigned int
forwards_coord3idx(const unsigned int idx,
                   const CartesianMaps &gbxmaps,
                   Superdrop &drop);

KOKKOS_FUNCTION unsigned int
backwards_coord1idx(const unsigned int idx,
                    const CartesianMaps &gbxmaps,
                    Superdrop &superdrop);
KOKKOS_FUNCTION unsigned int
forwards_coord1idx(const unsigned int idx,
                   const CartesianMaps &gbxmaps,
                   Superdrop &drop);

KOKKOS_FUNCTION unsigned int
backwards_coord2idx(const unsigned int idx,
                    const CartesianMaps &gbxmaps,
                    Superdrop &superdrop);
KOKKOS_FUNCTION unsigned int
forwards_coord2idx(const unsigned int idx,
                   const CartesianMaps &gbxmaps,
                   Superdrop &drop);

KOKKOS_FUNCTION void
check_inbounds_or_outdomain(const unsigned int idx,
                            const Kokkos::pair<double, double> bounds,
                            const double coord)
/* raise error if superdrop not either out of domain 
or within bounds (ie. lower_bound <= coord < upper_bound) */
{
  const bool bad_gbxindex((idx != outofbounds_gbxindex()) &&
                          ((coord < bounds.first) | (coord >= bounds.second)));
                          
  assert((!bad_gbxindex) && "SD not in previous gbx nor a neighbour."
                            " Try reducing the motion timestep to"
                            " satisfy CFL criteria, or use "
                            " 'update_ifoutside' to update sd_gbxindex");
}

int flag_sdgbxindex(const unsigned int idx,
                    const Kokkos::pair<double, double> bounds,
                    const double coord)
/* returns flag to keep idx the same (flag = 0) or
update to forwards (flag = 1) or backwards (flag = 2)
neighbour. Flag = 0 if idx is out of domain value or 
if coord lies within bounds = {lowerbound, upperbound}.
(Note: lower bound inclusive and upper bound exclusive,
ie. lowerbound <= coord < upperbound).
Flag = 1 if coord < lowerbound, indicating idx should 
be updated to backwards neighbour.
Flag = 2 if coord >= upperbound, indicating idx should 
be updated to forwards neighbour. */
{
  if (idx == outofbounds_gbxindex())
  {
    return 0; // maintian idx that is already out of domain
  }
  else if (coord < bounds.first) // lowerbound
  {
    return 1; // idx -> backwards_neighbour
  }
  else if (coord >= bounds.second) // upperbound
  {
    return 2; // idx -> forwards_neighbour
  }
  else
  {
    return 0; // maintain idx if coord within bounds
  }
}

KOKKOS_FUNCTION unsigned int
update_if_coord3nghbr(const CartesianMaps &gbxmaps,
                        unsigned int idx,
                        Superdrop &drop)
/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord3 direction. 
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
Note: backwards/forwards functions may change the 
superdroplet's attributes e.g. if it leaves the domain. */
{
  const int flag(flag_sdgbxindex(idx, gbxmaps.coord3bounds(idx),
                                 drop.get_coord3())); // if value != 0 idx needs to change
  switch (flag)
  {
  case 1:
    idx = backwards_coord3idx(idx, gbxmaps, drop);
    break;
  case 2:
    idx = forwards_coord3idx(idx, gbxmaps, drop);
    break;
  }
  return idx;
}

KOKKOS_FUNCTION unsigned int
update_if_coord1nghbr(const CartesianMaps &gbxmaps,
                          unsigned int idx,
                          Superdrop &drop)
/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord1 direction. 
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
Note: backwards/forwards functions may change the 
superdroplet's attributes e.g. if it leaves the domain. */
{
  const int flag(flag_sdgbxindex(idx, gbxmaps.coord1bounds(idx),
                                 drop.get_coord1())); // if value != 0 idx needs to change
  switch (flag)
  {
  case 1:
    idx = backwards_coord1idx(idx, gbxmaps, drop);
    break;
  case 2:
     idx = forwards_coord1idx(idx, gbxmaps, drop);
    break;
  }
  return idx;
}

KOKKOS_FUNCTION unsigned int
update_if_coord2nghbr(const CartesianMaps &gbxmaps,
                          unsigned int idx,
                          Superdrop &drop)
/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord2 direction. 
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
Note: backwards/forwards functions may change the 
superdroplet's attributes e.g. if it leaves the domain. */
{
  const int flag(flag_sdgbxindex(idx, gbxmaps.coord2bounds(idx),
                                 drop.get_coord2())); // if value != 0 idx needs to change
  switch (flag)
  {
  case 1:
    idx = backwards_coord2idx(idx, gbxmaps, drop);
    break;
  case 2:
     idx = forwards_coord2idx(idx, gbxmaps, drop);
    break;
  }
  return idx;
}

KOKKOS_FUNCTION unsigned int
backwards_coord3idx(const unsigned int idx,
                 const CartesianMaps &gbxmaps,
                 Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox
in backwards coord3 (z) direction and to update superdrop
coord3 if superdrop has exceeded the z lower domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord3backward(idx));

  const unsigned int incre(1); // increment
  if (at_cartesiandomainboundary(idx, incre, gbxmaps.get_ndim(0))) // drop was at lower z edge of domain (now moving below it)
  {
    const double lim1 = gbxmaps.coord3bounds(nghbr).second; // upper lim of backward neighbour
    const double lim2 = gbxmaps.coord3bounds(idx).first;    // lower lim of current gbx
    drop.set_coord3(coord3_beyondz(drop.get_coord3(), lim1, lim2));
  }

  return nghbr; // gbxindex of z backwards (down) neighbour
};

KOKKOS_FUNCTION unsigned int
forwards_coord3idx(const unsigned int idx,
                 const CartesianMaps &gbxmaps,
                 Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox in
forwards coord3 (z) direction and to update superdrop coord3
if superdrop has exceeded the z upper domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord3forward(idx));

  const unsigned int incre(1);                                             // increment
  if (at_cartesiandomainboundary(idx + incre, incre, gbxmaps.get_ndim(0))) // drop was upper z edge of domain (now moving above it)
  {
    const double lim1 = gbxmaps.coord3bounds(nghbr).first; // lower lim of forward neighbour
    const double lim2 = gbxmaps.coord3bounds(idx).second;  // upper lim of current gbx
    drop.set_coord3(coord3_beyondz(drop.get_coord3(), lim1, lim2));
  }

  return nghbr; // gbxindex of z forwards (up) neighbour
};

KOKKOS_FUNCTION unsigned int
backwards_coord1idx(const unsigned int idx,
                    const CartesianMaps &gbxmaps,
                    Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox
in backwards coord1 (x) direction and to update superdrop
coord1 if superdrop has exceeded the x back domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord1backward(idx));

  const auto ndims(gbxmaps.get_ndims());
  const unsigned int incre(ndims(0));                   // increment
  if (at_cartesiandomainboundary(idx, incre, ndims(1))) // at lower x edge of domain
  {
    const double lim1 = gbxmaps.coord1bounds(nghbr).second; // upper lim of backward neigghbour
    const double lim2 = gbxmaps.coord1bounds(idx).first;    // lower lim of current gbx
    drop.set_coord1(coord1_beyondx(drop.get_coord1(), lim1, lim2));
  }

  return nghbr; // gbxindex of x backwards (behind) neighbour
};

KOKKOS_FUNCTION unsigned int
forwards_coord1idx(const unsigned int idx,
                    const CartesianMaps &gbxmaps,
                    Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox
in forwards coord1 (x) direction and to update superdrop
coord1 if superdrop has exceeded the x front domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord1forward(idx));

  const auto ndims(gbxmaps.get_ndims());
  const unsigned int incre(ndims(0));                           // increment
  if (at_cartesiandomainboundary(idx + incre, incre, ndims(1))) // at lower x edge of domain
  {
    const double lim1 = gbxmaps.coord1bounds(nghbr).first; // lower lim of forward nghbour
    const double lim2 = gbxmaps.coord1bounds(idx).second;  // upper lim of gbx
    drop.set_coord1(coord1_beyondx(drop.get_coord1(), lim1, lim2));
  }

  return nghbr; // gbxindex of x forwards (infront) neighbour
};

KOKKOS_FUNCTION unsigned int
backwards_coord2idx(const unsigned int idx,
                    const CartesianMaps &gbxmaps,
                    Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox
in backwards coord2 (y) direction and to update superdrop
coord2 if superdrop has exceeded the y leftmost domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord2backward(idx));

  const auto ndims(gbxmaps.get_ndims());
  const unsigned int incre(ndims(0) * ndims(1));        // no. gridboxes in z direction * no. gridboxes in x direction
  if (at_cartesiandomainboundary(idx, incre, ndims(2))) // at lower y edge of domain
  {
    const double lim1 = gbxmaps.coord2bounds(nghbr).second; // upper lim of backward nghbour
    const double lim2 = gbxmaps.coord2bounds(idx).first;    // lower lim of gbx
    drop.set_coord2(coord2_beyondy(drop.get_coord2(), lim1, lim2));
  }

  return nghbr; // gbxindex of y backwards (left) neighbour
};

KOKKOS_FUNCTION unsigned int
forwards_coord2idx(const unsigned int idx,
                   const CartesianMaps &gbxmaps,
                   Superdrop &drop)
/* function to return gbxindex of neighbouring gridbox
in forwards coord2 (y) direction and to update superdrop
coord2 if superdrop has exceeded the y rightmost domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord2forward(idx));
  
  const auto ndims(gbxmaps.get_ndims());
  const unsigned int incre(ndims(0) * ndims(1));                            // no. gridboxes in z direction * no. gridboxes in x direction
  if (at_cartesiandomainboundary(idx + incre, incre, ndims(2)))             // at upper y edge of domain
  {
    const double lim1 = gbxmaps.coord2bounds(nghbr).first; // lower lim of forward nghbour
    const double lim2 = gbxmaps.coord2bounds(idx).second;  // upper lim of gbx
    drop.set_coord2(coord2_beyondy(drop.get_coord2(), lim1, lim2));
  }

  return nghbr; // gbxindex of y forwards (right) neighbour
};

#endif // CARTESIANMOTION_HPP