/*
 * ----- CLEO -----
 * File: cartesianmotion_withreset.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 19th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 19th December 2023
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
 * cartesian domain with finite/periodic boundary
 * conditions and reset of superdroplets that leave
 * the domain through the lower domain boundary
 */

#ifndef CARTESIANMOTION_WITHRESET_HPP
#define CARTESIANMOTION_WITHRESET_HPP

#include <functional>
#include <cassert>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "./cartesianmotion.hpp"
#include "./cartesianmaps.hpp"
#include "./cartesianboundaryconds.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"
#include "gridboxes/predcorr.hpp"

KOKKOS_FUNCTION unsigned int
change_to_coord3nghb_withreset(const CartesianMaps &gbxmaps,
                      unsigned int idx,
                      Superdrop &drop);

template <VelocityFormula TV>
struct CartesianMotionWithReset
/* satisfies motion concept for motion of a superdroplet
using a predictor-corrector method to update a superdroplet's
coordinates and then updating it's sdgbxindex using the
UpdateSdgbxindex struct for a cartesian domain */
{
  const unsigned int interval; // integer timestep for movement
  PredCorrMotion<CartesianMaps, TV> superdrop_coords;

  CartesianMotionWithReset(const unsigned int motionstep,
                           const std::function<double(unsigned int)> int2time,
                           const TV i_terminalv)
      : interval(motionstep),
        superdrop_coords(interval, int2time, i_terminalv) {}

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
  superdrop_gbx(const unsigned int gbxindex,
                const CartesianMaps &gbxmaps,
                Superdrop &drop) const
  /* function satisfies requirements of
  "superdrop_gbx" in the motion concept to update a
  superdroplet if it should move between gridboxes in a
  cartesian domain. For each direction (z, then x, then y),
  superdrop coord is compared to gridbox bounds given by gbxmaps
  for the current gbxindex 'idx'. If superdrop coord lies outside
  bounds, forward or backward neighbour functions are called as
  appropriate to update sdgbxindex (and possibly other attributes) */
  {
    unsigned int idx(gbxindex);

    idx = change_to_coord3nghbr_withreset(gbxmaps, idx, drop);
    check_inbounds_or_outdomain(idx, gbxmaps.coord3bounds(idx),
                                drop.get_coord3());

    idx = change_to_coord1nghbr(gbxmaps, idx, drop);
    check_inbounds_or_outdomain(idx, gbxmaps.coord1bounds(idx),
                                drop.get_coord1());

    idx = change_to_coord2nghbr(gbxmaps, idx, drop);
    check_inbounds_or_outdomain(idx, gbxmaps.coord2bounds(idx),
                                drop.get_coord2());

    drop.set_sdgbxindex(idx);
  }
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION unsigned int
change_to_coord3nghbr_withreset(const CartesianMaps &gbxmaps,
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
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord3bounds(idx),
                                    drop.get_coord3()); // if value != 0 idx needs to change
  switch (flag)
  {
  case 1:
    idx = change_to_backwards_coord3nghbr_withreset(idx, gbxmaps, drop);
    break;
  case 2:
    idx = change_to_forwards_coord3nghbr(idx, gbxmaps, drop);
    break;
  }
  return idx;
}

#endif // CARTESIANMOTION_WITHRESET_HPP

