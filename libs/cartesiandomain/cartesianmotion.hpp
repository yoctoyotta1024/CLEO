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
#include "gridboxes/predictorcorrector.hpp"

struct UpdateSdgbxindex
/* struct contanies operator to satisfiy
requirements of "update_superdrop_gbxindex"
in the motion concept. Operator updates
superdroplet sdgbxindex in a cartesian domain */
{
  KOKKOS_INLINE_FUNCTION void
  operator()(const unsigned int gbxindex,
             const CartesianMaps &gbxmaps,
             Superdrop &drop) const {}
};

template <VelocityFormula TerminalVelocity>
struct CartesianMotion
/* satisfies motion concept for motion of a superdroplet
using a predictor-corrector method to update a superdroplet's
coordinates and then updating it's sdgbxindex using the
UpdateSdgbxindex struct for a cartesian domain */
{
  const unsigned int interval; // integer timestep for movement
  
  PredCorrMotion<CartesianMaps, TerminalVelocity>
      update_superdrop_coords;

  UpdateSdgbxindex update_superdrop_gbxindex;

  CartesianMotion(const unsigned int motionstep,
                  const std::function<double(int)> int2time,
                  const TerminalVelocity i_terminalv)
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
};

#endif // CARTESIANMOTION_HPP