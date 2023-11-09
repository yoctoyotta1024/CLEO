/*
 * ----- CLEO -----
 * File: cartesianmotion.hpp
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
#include "./boundaryconditions.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"
#include "gridboxes/predcorr.hpp"

KOKKOS_FUNCTION void
cartesian_update_superdrop_gbxindex(const unsigned int gbxindex,
                                    const CartesianMaps &gbxmaps,
                                    Superdrop &drop);

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
                  const std::function<double(int)> int2time,
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
  "update_superdrop_gbxindex" in the motion concept.
  calls function for updating superdroplet sdgbxindex
  in a cartesian domain */
  {
    cartesian_update_superdrop_gbxindex(gbxindex, gbxmaps, drop);
  }
};

#endif // CARTESIANMOTION_HPP