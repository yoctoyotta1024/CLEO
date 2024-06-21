/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: gridboxes
 * Created Date: Friday 3rd May 2024
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


/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: gridboxes
 * Created Date: Tuesday 19th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 30th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * -----
 * File Description:
 * Generic struct satisfying Motion concept for
 * a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * updating gbx according to templated functions
 */

#ifndef LIBS_GRIDBOXES_PREDCORRMOTION_HPP_
#define LIBS_GRIDBOXES_PREDCORRMOTION_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <functional>

#include "gridboxes/predcorr.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

/* satisfies motion concept for motion of a superdroplet
using a predictor-corrector method to update a superdroplet's
coordinates and then updating it's sdgbxindex using
the appropriate templated type */
template <GridboxMaps GbxMaps, VelocityFormula TV, typename ChangeToNghbr, typename CheckBounds>
struct PredCorrMotion {
  const unsigned int interval;  // integer timestep for movement
  PredCorr<GbxMaps, TV> superdrop_coords;
  ChangeToNghbr change_if_nghbr;
  CheckBounds check_bounds;

  PredCorrMotion(const unsigned int motionstep, const std::function<double(unsigned int)> int2time,
                 const TV i_terminalv, ChangeToNghbr i_change_if_nghbr, CheckBounds i_check_bounds)
      : interval(motionstep),
        superdrop_coords(interval, int2time, i_terminalv),
        change_if_nghbr(i_change_if_nghbr),
        check_bounds(i_check_bounds) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const {
    return ((t_sdm / interval) + 1) * interval;
  }

  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int t_sdm) const { return t_sdm % interval == 0; }

  /* function satisfies requirements of
  "superdrop_gbx" in the motion concept to update a
  superdroplet if it should move between gridboxes.
  For each direction (coord3, then coord1, then coord2),
  superdrop and idx may be changed if superdrop coord
  lies outside bounds of gridbox in that direction */
  KOKKOS_INLINE_FUNCTION void superdrop_gbx(const unsigned int gbxindex,
                                            const CartesianMaps &gbxmaps, Superdrop &drop) const {
    auto idx = (unsigned int)change_if_nghbr.coord3(gbxmaps, gbxindex, drop);
    check_bounds(idx, gbxmaps.coord3bounds(idx), drop.get_coord3());

    idx = change_if_nghbr.coord1(gbxmaps, idx, drop);
    check_bounds(idx, gbxmaps.coord1bounds(idx), drop.get_coord1());

    idx = change_if_nghbr.coord2(gbxmaps, idx, drop);
    check_bounds(idx, gbxmaps.coord2bounds(idx), drop.get_coord2());

    assert((drop.get_sdgbxindex() == idx) && "sdgbxindex not concordant with supposed idx");
  }
};

#endif  // LIBS_GRIDBOXES_PREDCORRMOTION_HPP_
