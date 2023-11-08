/*
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: gridboxes
 * Created Date: Monday 16th October 2023
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
 * Change in a superdroplet's coords using predictor-
 * corrector method for motion of superdroplet given
 * a formula for its terminal velocity and the wind
 * velocity obtained via a simple linear interpolation.
 * Methods follows equations in Grabowski et al. 2018
 */

#ifndef PREDCORRMOTION_HPP
#define PREDCORRMOTION_HPP

#include <cassert>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "./cfl_criteria.hpp"
#include "./gridboxmaps.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/state.hpp"
#include "superdrops/terminalvelocity.hpp"

KOKKOS_FUNCTION
double interpolation(const Kokkos::pair<double, double> bounds,
                     const Kokkos::pair<double, double> vel,
                     const double sdcoord);
/* Given [X = z,x or y] wind velocity component, vel, that is
defined on the faces of a gridbox at {lower, upper} [X] bounds,
return wind at [X] coord. Method is 'simple' linear interpolation
from Grabowski et al. (2018). coord use in interpolation is
limited to lower_bound <= coord <= upper_bound. */

template <GridboxMaps GbxMaps, VelocityFormula TerminalVelocity>
struct PredCorrMotion
/* change in coordinates calculated by predictor
corrector method with the wind velocity obtained via
a simple linear interpolation. Methods follows
equations in Grabowski et al. 2018 */
{
private:
  const unsigned int interval; // integer timestep for movement
  const double delt;           // equivalent of interval as dimensionless time
  const TerminalVelocity terminalv; // returns terminal velocity given a superdroplet

  KOKKOS_INLINE_FUNCTION
  double interp_wvel(const unsigned int gbxindex,
                     const GbxMaps &gbxmaps,
                     const State &state,
                     const double coord3) const
  /* method to interpolate coord3 wind velocity component (w)
  defined on coord3 faces of a gridbox to a superdroplet's
  coordinates at (coord3, coord1, coord2) */
  {
    return interpolation(gbxmaps.coord3bounds(gbxindex),
                         state.wvel, coord3);
  }

  KOKKOS_INLINE_FUNCTION
  double interp_uvel(const unsigned int gbxindex,
                     const GbxMaps &gbxmaps,
                     const State &state,
                     const double coord1) const
  /* method to interpolate coord1 wind velocity component (u)
  defined on coord1 faces of a gridbox to a superdroplet's
  coordinates at (coord3, coord1, coord2) */
  {
    return interpolation(gbxmaps.coord1bounds(gbxindex),
                         state.uvel, coord1);
  }

  KOKKOS_INLINE_FUNCTION
  double interp_vvel(const unsigned int gbxindex,
                     const GbxMaps &gbxmaps,
                     const State &state,
                     const double coord2) const
  /* method to interpolate coord2 wind velocity component (v)
  defined on coord2 faces of a gridbox to a superdroplet's
  coordinates at (coord3, coord1, coord2) */
  {
    return interpolation(gbxmaps.coord2bounds(gbxindex),
                         state.vvel, coord2);
  }

  KOKKOS_FUNCTION
  double delta_coord3(const unsigned int gbxindex,
                      const GbxMaps &gbxmaps,
                      const State &state,
                      const Superdrop &drop) const
  {
    double coord3(drop.get_coord3());

    const double terminal = terminalv(drop);

    /* corrector velocities based on predicted coords */
    double vel3 = interp_wvel(gbxindex, gbxmaps, state, coord3);
    vel3 -= terminal;

    /* predictor coords given velocity at previous coords */
    coord3 += vel3 * delt; // move by w wind + terminal velocity

    /* corrector velocities based on predicted coords */
    double corrvel3 = interp_wvel(gbxindex, gbxmaps, state, coord3);
    corrvel3 -= terminal;

    /* predicted-corrected change to superdrop coords */
    const double delta3((vel3 + corrvel3) * (delt / 2));

    return delta3;
  }

  KOKKOS_FUNCTION
  double delta_coord1(const unsigned int gbxindex,
                      const GbxMaps &gbxmaps,
                      const State &state,
                      const Superdrop &drop) const
  {
    double coord1(drop.get_coord1());

    /* corrector velocities based on predicted coords */
    const double vel1 = interp_uvel(gbxindex, gbxmaps, state, coord1);

    /* predictor coords given velocity at previous coords */
    coord1 += vel1 * delt; // move by u wind

    /* corrector velocities based on predicted coords */
    const double corrvel1 = interp_uvel(gbxindex, gbxmaps, state, coord1);

    /* predicted-corrected change to superdrop coords */
    const double delta1((vel1 + corrvel1) * (delt / 2));

    return delta1;
  }

  KOKKOS_FUNCTION
  double delta_coord2(const unsigned int gbxindex,
                      const GbxMaps &gbxmaps,
                      const State &state,
                      const Superdrop &drop) const
  {
    double coord2(drop.get_coord2());

    /* corrector velocities based on predicted coords */
    const double vel2 = interp_vvel(gbxindex, gbxmaps, state, coord2);

    /* predictor coords given velocity at previous coords */
    coord2 += vel2 * delt; // move by v wind

    /* corrector velocities based on predicted coords */
    const double corrvel2 = interp_vvel(gbxindex, gbxmaps, state, coord2);

    /* predicted-corrected change to superdrop coords */
    const double delta2((vel2 + corrvel2) * (delt / 2));

    return delta2;
  }

public:
  PredCorrMotion(const unsigned int motionstep,
                 const std::function<double(int)> int2time)
      : interval(motionstep),
        delt(int2time(interval)),
        terminalv(TerminalVelocity{}) {}

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

  KOKKOS_FUNCTION
  void update_superdrop_coords(const unsigned int gbxindex,
                               const GbxMaps &gbxmaps,
                               const State &state,
                               Superdrop &drop) const
  /* Uses predictor-corrector method to forward timestep
  a superdroplet's coordinates using the interpolated
  wind velocity from a gridbox's state */
  {
    /* Use predictor-corrector method to get change in SD coords */
    const double delta3 = delta_coord3(gbxindex, gbxmaps, state, drop);
    const double delta1 = delta_coord1(gbxindex, gbxmaps, state, drop);
    const double delta2 = delta_coord2(gbxindex, gbxmaps, state, drop);

    /* CFL check on predicted change to SD coords */
    cfl_criteria(gbxmaps, gbxindex, delta3, delta1, delta2);

    /* update SD coords */
    drop.increment_coords(delta3, delta1, delta2);
  }
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION
double interpolation(const Kokkos::pair<double, double> bounds,
                     const Kokkos::pair<double, double> vel,
                     const double sdcoord)
/* Given [X = z,x or y] wind velocity component, vel, that is
defined on the faces of a gridbox at {lower, upper} [X] bounds,
return wind at [X] coord. Method is 'simple' linear interpolation
from Grabowski et al. (2018). coord use in interpolation is
limited to lower_bound <= coord <= upper_bound. */
{
  const double coord(Kokkos::fmin(bounds.second,
                                  Kokkos::fmax(bounds.first,
                                               sdcoord))); // limit coord to within bounds

  const double alpha((coord - bounds.first) /
                     (bounds.second - bounds.first));

  const double interp(alpha * vel.second +
                      (1 - alpha) * vel.first); // simple linear interpolation

  return interp;
}

#endif // PREDCORRMOTION_HPP
