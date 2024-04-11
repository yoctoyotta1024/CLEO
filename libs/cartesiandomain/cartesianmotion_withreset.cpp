/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesianmotion_withreset.cpp
 * Project: cartesiandomain
 * Created Date: Tuesday 19th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 11th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Some functionality for motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * the sdgbxindex updated accordingly for a
 * cartesian domain with finite/periodic boundary
 * conditions and reset of superdroplets that leave
 * the domain through coord3 domain boundaries
 */

#include "./cartesianmotion_withreset.hpp"

KOKKOS_FUNCTION unsigned int change_to_backwards_coord3nghbr_withreset(
    const ResetSuperdrop &reset_superdrop, const unsigned int idx, const CartesianMaps &gbxmaps,
    Superdrop &superdrop);

KOKKOS_FUNCTION unsigned int change_to_forwards_coord3nghbr_withreset(
    const ResetSuperdrop &reset_superdrop, const unsigned int idx, const CartesianMaps &gbxmaps,
    Superdrop &superdrop);

/* return updated value of gbxindex in case superdrop should
move to neighbouring gridbox in coord3 direction.
Funciton changes value of idx if flag != 0,
if flag = 1 idx updated to backwards neighbour gbxindex.
if flag = 2 idx updated to forwards neighbour gbxindex.
_Note:_ backwards/forwards functions may change the
superdroplet's attributes e.g. if it leaves the domain. */
KOKKOS_FUNCTION unsigned int change_if_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                                             const CartesianMaps &gbxmaps,
                                                             unsigned int idx, Superdrop &drop) {
  const auto flag = flag_sdgbxindex(idx, gbxmaps.coord3bounds(idx),
                                    drop.get_coord3());  // if value != 0 idx needs to change
  switch (flag) {
    case 1:
      idx = change_to_backwards_coord3nghbr_withreset(reset_superdrop, idx, gbxmaps, drop);
      break;
    case 2:
      idx = change_to_forwards_coord3nghbr_withreset(reset_superdrop, idx, gbxmaps, drop);
      break;
  }
  return idx;
}

/* function to return gbxindex of neighbouring gridbox
in backwards coord3 (z) direction and to update superdrop
if its coord3 has exceeded the z lower domain boundary */
KOKKOS_FUNCTION unsigned int change_to_backwards_coord3nghbr_withreset(
    const ResetSuperdrop &reset_superdrop, const unsigned int idx, const CartesianMaps &gbxmaps,
    Superdrop &drop) {
  auto nghbr = (unsigned int)gbxmaps.coord3backward(idx);

  const auto incre = (unsigned int)1;  // increment
  // drop was at lower z edge of domain (now moving below it)
  if (beyond_domainboundary(idx, incre, gbxmaps.get_ndim(0))) {
    nghbr = reset_superdrop(gbxmaps, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of z backwards (down) neighbour
}

/* function to return gbxindex of neighbouring gridbox in
forwards coord3 (z) direction and to update superdrop coord3
if superdrop has exceeded the z upper domain boundary */
KOKKOS_FUNCTION unsigned int change_to_forwards_coord3nghbr_withreset(
    const ResetSuperdrop &reset_superdrop, const unsigned int idx, const CartesianMaps &gbxmaps,
    Superdrop &drop) {
  auto nghbr = (unsigned int)gbxmaps.coord3forward(idx);

  const auto incre = (unsigned int)1;  // increment
  // drop was upper z edge of domain (now moving above it)
  if (beyond_domainboundary(idx + incre, incre, gbxmaps.get_ndim(0))) {
    nghbr = reset_superdrop(gbxmaps, drop);
  }

  drop.set_sdgbxindex(nghbr);
  return nghbr;  // gbxindex of z forwards (up) neighbour
}
