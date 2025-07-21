/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: predcorr.hpp
 * Project: gridboxes
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Change in a superdroplet's coords using predictor-
 * corrector method for motion of superdroplet given
 * a formula for its terminal velocity and the wind
 * velocity obtained via a simple linear interpolation.
 * Methods follows equations in Grabowski et al. 2018
 */

#ifndef LIBS_GRIDBOXES_PREDCORR_HPP_
#define LIBS_GRIDBOXES_PREDCORR_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <cassert>
#include <functional>

#include "../kokkosaliases.hpp"
#include "gridboxes/cfl_criteria.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

/* Given [X = z,x or y] wind velocity component, vel, that is
defined on the faces of a gridbox at {lower, upper} [X] bounds,
return wind at [X] coord. Method is 'simple' linear interpolation
from Grabowski et al. (2018). coord use in interpolation is
limited to lower_bound <= coord <= upper_bound. */
KOKKOS_FUNCTION
double interpolation(const Kokkos::pair<double, double> bounds,
                     const Kokkos::pair<double, double> vel, const double sdcoord);

/* change in coordinates calculated by predictor
corrector method with the wind velocity obtained via
a simple linear interpolation. Methods follows
equations in Grabowski et al. 2018 */
template <GridboxMaps GbxMaps, VelocityFormula TV>
struct PredCorr {
 private:
  const double delt;   // equivalent of motionstep as dimensionless time
  const TV terminalv;  // returns terminal velocity given a superdroplet

  /* method to interpolate coord3 wind velocity component (w)
  defined on coord3 faces of a gridbox to a superdroplet's
  coordinates at (coord3, coord1, coord2) */
  KOKKOS_INLINE_FUNCTION
  double interp_wvel(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                     const double coord3) const {
    return interpolation(gbxmaps.coord3bounds(gbxindex), state.wvel, coord3);
  }

  /* method to interpolate coord1 wind velocity component (u)
  defined on coord1 faces of a gridbox to a superdroplet's
  coordinates at (coord3, coord1, coord2) */
  KOKKOS_INLINE_FUNCTION
  double interp_uvel(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                     const double coord1) const {
    return interpolation(gbxmaps.coord1bounds(gbxindex), state.uvel, coord1);
  }

  /* method to interpolate coord2 wind velocity component (v)
  defined on coord2 faces of a gridbox to a superdroplet's
  coordinates at (coord3, coord1, coord2) */
  KOKKOS_FUNCTION
  double interp_vvel(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                     const double coord2) const {
    return interpolation(gbxmaps.coord2bounds(gbxindex), state.vvel, coord2);
  }

  KOKKOS_FUNCTION
  double delta_coord3(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                      const Superdrop &drop) const {
    auto coord3 = drop.get_coord3();

    const auto terminal = terminalv(drop);

    /* corrector velocities based on predicted coords */
    auto vel3 = interp_wvel(gbxindex, gbxmaps, state, coord3);
    vel3 -= terminal;

    /* predictor coords given velocity at previous coords */
    coord3 += vel3 * delt;  // move by w wind + terminal velocity

    /* corrector velocities based on predicted coords */
    auto corrvel3 = interp_wvel(gbxindex, gbxmaps, state, coord3);
    corrvel3 -= terminal;

    /* predicted-corrected change to superdrop coords */
    const auto delta3 = double{(vel3 + corrvel3) * (delt / 2)};

    return delta3;
  }

  KOKKOS_FUNCTION
  double delta_coord1(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                      const Superdrop &drop) const {
    auto coord1 = drop.get_coord1();

    /* corrector velocities based on predicted coords */
    const auto vel1 = interp_uvel(gbxindex, gbxmaps, state, coord1);

    /* predictor coords given velocity at previous coords */
    coord1 += vel1 * delt;  // move by u wind

    /* corrector velocities based on predicted coords */
    const auto corrvel1 = interp_uvel(gbxindex, gbxmaps, state, coord1);

    /* predicted-corrected change to superdrop coords */
    const auto delta1 = double{(vel1 + corrvel1) * (delt / 2)};

    return delta1;
  }

  KOKKOS_FUNCTION
  double delta_coord2(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                      const Superdrop &drop) const {
    auto coord2 = drop.get_coord2();

    /* corrector velocities based on predicted coords */
    const auto vel2 = interp_vvel(gbxindex, gbxmaps, state, coord2);

    /* predictor coords given velocity at previous coords */
    coord2 += vel2 * delt;  // move by v wind

    /* corrector velocities based on predicted coords */
    const auto corrvel2 = interp_vvel(gbxindex, gbxmaps, state, coord2);

    /* predicted-corrected change to superdrop coords */
    const auto delta2 = double{(vel2 + corrvel2) * (delt / 2)};

    return delta2;
  }

 public:
  PredCorr(const unsigned int motionstep, const std::function<double(unsigned int)> int2time,
           const TV i_terminalv)
      : delt(int2time(motionstep)), terminalv(i_terminalv) {}

  /* operator for use in the "superdrop_coords" function of the PredCorrMotion struct.
  Operator uses predictor-corrector method to return the change in
  a superdroplet's coordinates from a forward timestep of motion using the
  interpolated wind velocity from a gridbox's state */
  KOKKOS_FUNCTION
  Superdrop operator()(const unsigned int gbxindex, const GbxMaps &gbxmaps, const State &state,
                       Superdrop &drop) const {
    /* Use predictor-corrector method to get change in SD coords */
    const auto delta3 = delta_coord3(gbxindex, gbxmaps, state, drop);
    const auto delta1 = delta_coord1(gbxindex, gbxmaps, state, drop);
    const auto delta2 = delta_coord2(gbxindex, gbxmaps, state, drop);

    /* CFL check on predicted change to SD coords */
    cfl_criteria(gbxmaps, gbxindex, delta3, delta1, delta2);

    /* update SD coords */
    drop.increment_coords(delta3, delta1, delta2);

    return drop;
  }
};

#endif  // LIBS_GRIDBOXES_PREDCORR_HPP_
