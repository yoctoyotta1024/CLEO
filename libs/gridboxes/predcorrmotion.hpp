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

#include "../cleoconstants.hpp"
#include "gridboxes/predcorr.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

/*
satisfies motion concept for motion of a superdroplet using a predictor-corrector method with
a constant timestep ("interval") to update a superdroplet's coordinates and then updating it's
sdgbxindex using the appropriate templated type

Special case: If timestep interval is largest possible unsigned integer, on_step never returns true.
*/
template <GridboxMaps GbxMaps, VelocityFormula TV, typename CheckBounds>
struct PredCorrMotion {
  const unsigned int interval;  // integer timestep for movement
  PredCorr<GbxMaps, TV> predcorr;
  CheckBounds check_bounds;

  PredCorrMotion(const unsigned int motionstep, const std::function<double(unsigned int)> int2time,
                 const TV i_terminalv, CheckBounds i_check_bounds)
      : interval(motionstep),
        predcorr(interval, int2time, i_terminalv),
        check_bounds(i_check_bounds) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const {
    return ((t_sdm / interval) + 1) * interval;
  }

  /**
   * @brief Returns true if motion should perform an on-step action.
   *
   * Special case: If interval is largest possible unsigned integer, on_step never returns true.
   *
   * @param t_sdm The current time step.
   * @return True if the current time step is a multiple of the interval.
   */
  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int t_sdm) const {
    return (t_sdm % interval == 0) && (interval != LIMITVALUES::uintmax);
  }

  /* function satisfies requirements of the "superdrop_coords" function in the motion
  concept. Operator uses predictor-corrector method to obtain the change in the coordinates
  from a forward timestep of the motion and then increments the superdroplet coordinates
  accordingly */
  KOKKOS_INLINE_FUNCTION
  void superdrop_coords(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                        Superdrop &drop) const {
    /* change in SD coords: (coord3, coord1, coord2) */
    drop = predcorr(gbxindex, gbxmaps, state, drop);
  }

  /* function satisfies requirements of "superdrop_gbx" in the motion
  concept to update a superdroplet if it should move between gridboxes
  (or out of domain). Function also called check_bounds to check
  superdroplet is indeed in correct gridbox after update. */
  KOKKOS_INLINE_FUNCTION void superdrop_gbx(const unsigned int gbxindex, const GbxMaps &gbxmaps,
                                            Superdrop &drop) const {
    auto coord3 = drop.get_coord3();
    auto coord1 = drop.get_coord1();
    auto coord2 = drop.get_coord2();
    const auto idx = gbxmaps.get_local_bounding_gridbox_index(gbxindex, coord3, coord1,
                                                        coord2);  // drop_coords may change(!)

    // Sets the updated superdroplet coordinates and gridbox index
    drop.set_coords(coord3, coord1, coord2);
    drop.set_sdgbxindex(idx);

    // If the index is non-local return
    // For superdrops going to other processes checks will be perfomed in the receiver
    // For out of bounds index nothing will be done
    if (idx >= gbxmaps.get_local_ngridboxes()) return;

    // Checks that the drop coordinates match the ones of the gridbox
    check_bounds(idx, gbxmaps.coord3bounds(idx), drop.get_coord3());
    check_bounds(idx, gbxmaps.coord1bounds(idx), drop.get_coord1());
    check_bounds(idx, gbxmaps.coord2bounds(idx), drop.get_coord2());
  }
};

#endif  // LIBS_GRIDBOXES_PREDCORRMOTION_HPP_
