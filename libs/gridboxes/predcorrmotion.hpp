/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: gridboxes
 * Created Date: Tuesday 19th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 9th July 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
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
    auto drop_coords =
        std::array<double, 3>{drop.get_coord3(), drop.get_coord1(), drop.get_coord2()};
    const auto domain_decomposition = gbxmaps.get_domain_decomposition();
    const auto idx = domain_decomposition.get_local_bounding_gridbox(drop_coords);

    // Sets the updated superdroplet coordinates and gridbox index
    drop.set_coord3(drop_coords[0]);
    drop.set_coord1(drop_coords[1]);
    drop.set_coord2(drop_coords[2]);
    drop.set_sdgbxindex(idx);

    // If the index is non-local return
    // For superdrops going to other processes checks will be perfomed in the receiver
    // For out of bounds index nothing will be done
    if (idx >= domain_decomposition.get_total_local_gridboxes()) return;

    // Checks that the drop coordinates match the ones of the gridbox
    check_bounds(idx, gbxmaps.coord3bounds(idx), drop.get_coord3());
    check_bounds(idx, gbxmaps.coord1bounds(idx), drop.get_coord1());
    check_bounds(idx, gbxmaps.coord2bounds(idx), drop.get_coord2());

    assert((drop.get_sdgbxindex() == idx) && "sdgbxindex not concordant with supposed idx");
  }
};

#endif  // LIBS_GRIDBOXES_PREDCORRMOTION_HPP_
