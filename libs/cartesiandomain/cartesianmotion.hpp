/*
 * ----- CLEO -----
 * File: cartesianmotion.hpp
 * Project: cartesiandomain
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 8th November 2023
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

#include "./cartesianmaps.hpp"
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

KOKKOS_FUNCTION void
cartesian_update_superdrop_gbxindex(const unsigned int gbxindex,
                                    const CartesianMaps &gbxmaps,
                                    Superdrop &drop)
/* Updates superdroplet sdgbxindex in a cartesian domain.
For each direction (z, then x, then y), gbxmaps's forward and backward
get_neighbour functions are passed into update_superdrop_ifneighbour
along with superdroplet and the gridbox bounds for that direction.
(If coord not within bounds, update_superdrop_ifneighbour should call
appropriate get_neighbour function to update the superdroplet's
sd_gbxindex (and possibly other attributes if desired). After algorithm
for z, then x, then y directions are complete, resultant sd_gbxindex
is returned. */
{
  unsigned int current_gbxindex(gbxindex);

  current_gbxindex = update_ifneighbour(
      gbxmaps, zdown, zup,
      [](const Maps4GridBoxes &gbxmaps, const auto ii)
      { return gbxmaps.get_bounds_z(ii); },
      [](const Superdrop &superdrop)
      { return superdrop.coord3; },
      current_gbxindex, superdrop);

  current_gbxindex = update_ifneighbour(
      gbxmaps, xbehind, xinfront,
      [](const Maps4GridBoxes &gbxmaps,
         const unsigned int ii)
      { return gbxmaps.get_bounds_x(ii); },
      [](const Superdrop &superdrop)
      { return superdrop.coord1; },
      current_gbxindex, superdrop);

  current_gbxindex = update_ifneighbour(
      gbxmaps, yleft, yright,
      [](const Maps4GridBoxes &gbxmaps,
         const unsigned int ii)
      { return gbxmaps.get_bounds_y(ii); },
      [](const Superdrop &superdrop)
      { return superdrop.coord2; },
      current_gbxindex, superdrop);

  drop.set_sdgbxindex(current_gbxindex);
}

#endif // CARTESIANMOTION_HPP