/*
 * ----- CLEO -----
 * File: cartesianmotion.cpp
 * Project: cartesiandomain
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 9th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 */

#include "./cartesianmotion.hpp"

KOKKOS_FUNCTION void
check_inbounds_or_outdomain(const unsigned int idx,
                            const Kokkos::pair<double, double> bounds,
                            const double coord);
/* raise error if superdrop not out of domain or within bounds */

KOKKOS_FUNCTION int
flag_sdgbxindex(const unsigned int idx,
                const Kokkos::pair<double, double> bounds,
                const double coord);

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

KOKKOS_FUNCTION unsigned int
backwards_coord3(const unsigned int gbxindex,
                 const CartesianMaps &gbxmaps,
                 Superdrop &superdrop);

KOKKOS_FUNCTION unsigned int
forwards_coord3(const unsigned int gbxindex,
                const CartesianMaps &gbxmaps,
                Superdrop &drop);

KOKKOS_FUNCTION void
CartesianMotion::update_superdrop_gbxindex(const unsigned int gbxindex,
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

  idx = update_if_coord1nghbr(gbxmaps, idx, drop);

  idx = update_if_coord2nghbr(gbxmaps, idx, drop);

  // current_gbxindex = update_ifneighbour(
  //                   gbxmaps,
  //                   backwards_neighbour = xbehind,
  //                   forwards_neighbour = xinfront,
  //                   get_bounds = gbxmaps.get_bounds_x(ii),
  //                   get_sdcoord = superdrop.coord1,
  //                   idx,
  //                   superdrop);

  // current_gbxindex = update_ifneighbour(
  //                   gbxmaps,
  //                   backwards_neighbour =  yleft,
  //                   forwards_neighbour =yright,
  //                   get_bounds = gbxmaps.get_bounds_y(ii),
  //                   get_sdcoord = superdrop.coord2,
  //                   idx,
  //                   superdrop);

  drop.set_sdgbxindex(idx);
}

KOKKOS_FUNCTION void
check_inbounds_or_outdomain(const unsigned int idx,
                            const Kokkos::pair<double, double> bounds,
                            const double coord)
/* raise error if superdrop not out of domain or within bounds */
{
  const bool bad_gbxindex((idx != LIMITVALUES::uintmax) &&
                          ((coord < bounds.first) | (coord >= bounds.second)));

  assert((!bad_gbxindex) && "SD not in previous gbx nor a neighbour."
                            " Try reducing the motion timestep to"
                            " satisfy CFL criteria, or use "
                            " 'update_ifoutside' to update sd_gbxindex");
}

int flag_sdgbxindex(const unsigned int idx,
                    const Kokkos::pair<double, double> bounds,
                    const double coord)
/* Given bounds = {lowerbound, upperbound} of a gridbox with
index 'gbxindex', function determines if coord is within bounds
of that gridbox. (Note: lower bound inclusive, upper bound exclusive).
If coord not within bounds backwardsidx or forwardsidx function,
as appropriate, is used to return a neighbouring gridbox's index.
If coord lies within bounds, gbxindex is returned. If index is
already out of domain (ie. value is the maximum unsigned int),
return out of domain index */
{
  if (idx == LIMITVALUES::uintmax)
  {
    return 0; // ignore idx that is already out of domain
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
    return 0; // no change to idx if coord within bounds
  }
}

KOKKOS_FUNCTION unsigned int
update_if_coord3neighbour(const CartesianMaps &gbxmaps,
                        unsigned int idx,
                        Superdrop &drop)
/* For a given direction, pass {lower, upper} bounds into
update_superdrop_ifneighbour to get updated gbxindex and superdrop
(e.g. if superdroplet's coord fromsdcoord function lies outside of
bounds given gbxbounds using idx). Repeat until
superdroplet coord is within the bounds given by the idx,
or until superdrop leaves domain. */
{
  const int flag(flag_sdgbxindex(idx,
                                 gbxmaps.coord3bounds(idx),
                                 drop.get_coord3())); // value != 0 if sdgbxindex needs to change
  switch (flag)
  {
  case 1:
    idx = backwards_coord3(idx, gbxmaps, drop);
    break;
  case 2:
     idx = forwards_coord3(idx, gbxmaps, drop);
    break;
  }

  check_inbounds_or_outdomain(idx,
                              gbxmaps.coord3bounds(idx),
                              drop.get_coord3());

  return idx;
}

KOKKOS_FUNCTION unsigned int
backwards_coord3(const unsigned int idx,
                 const CartesianMaps &gbxmaps,
                 Superdrop &drop)
/* function to return sdgbxindex of neighbouring gridbox
in backwards coord3 (z) direction and to update superdrop
coord3 if superdrop has exceeded the z lower domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord3backward(idx));

  if (at_cartesiandomainboundary(idx, 1, gbxmaps.get_ndim(0))) // SD was at lower z edge of domain (now moving beyond it)
  {
    const double lim1 = gbxmaps.coord3bounds(nghbr).second;   // upper lim of backward nghbour
    const double lim2 = gbxmaps.coord3bounds(idx).first; // lower lim of gbx
    drop.set_coord3(coord3_beyondz(drop.get_coord3(), lim1, lim2));
  }

  return nghbr; // gbxindex of zdown_neighbour
};

KOKKOS_FUNCTION unsigned int
forwards_coord3(const unsigned int idx,
                 const CartesianMaps &gbxmaps,
                 Superdrop &drop)
/* function to return sdgbxindex of neighbouring gridbox in
forwards coord3 (z) direction and to update superdrop coord3
if superdrop has exceeded the z upper domain boundary */
{
  const unsigned int nghbr(gbxmaps.coord3forward(idx));

  if (at_cartesiandomainboundary(idx + 1, 1, gbxmaps.get_ndim(0))) // SD was upper z edge of domain (now moving above it)
  {
    const double lim1 = gbxmaps.coord3bounds(nghbr).first; // lower lim of forward nghbour
    const double lim2 = gbxmaps.coord3bounds(idx).second; // upper lim of gbx
    drop.set_coord3(coord3_beyondz(drop.get_coord3(), lim1, lim2));
  }

  return nghbr; // gbxindex of zup_neighbour
};