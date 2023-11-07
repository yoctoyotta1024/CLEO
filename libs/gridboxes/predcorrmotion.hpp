/*
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: gridboxes
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates given
 * a formula for its terminal velocity and the wind
 * velocity obtained via a simple linear interpolation.
 * Methods follows equations in Grabowski et al. 2018
 */


#ifndef PREDCORRMOTION_HPP
#define PREDCORRMOTION_HPP

#include <cassert>

#include <Kokkos_Core.hpp>

#include "./cfl_criteria.hpp"
#include "./gridboxmaps.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/state.hpp"

struct InterpolateWinds
/* method to interpolate (w, u, and v) wind velocity components
defined on (coord3, coord1 and coord2) faces of a gridbox
to a superdroplet's coordinates at (coord3, coord1, coord2) */
{
}

template <GridboxMaps GbxMaps>
struct PredCorrMotion
/* change in coordinates calculated by predictor
corrector method with the wind velocity obtained via 
a simple linear interpolation. Methods follows
equations in Grabowski et al. 2018 */
{
private:
  unsigned int interval;

  KOKKOS_INLINE_FUNCTION
  double delta_coord3(Superdrop &drop
                      InterpolateWinds &windinterp)
  {
    const double terminal = terminalv(drop);

    /* corrector velocities based on predicted coords */
    const double vel3 = windinterp.interp_wvel() - terminal;

    /* predictor coords given velocity at previous coords */
    windinterp.coord3 += vel3 * delt; // move by w wind + terminal velocity

    /* corrector velocities based on predicted coords */
    const double corrvel3 = windinterp.interp_wvel() - terminal;

    /* predicted-corrected change to superdrop coords */
    const double delta3((vel3 + corrvel3) * (delt / 2));
  }

  KOKKOS_INLINE_FUNCTION
  double delta_coord1(Superdrop &drop
                      InterpolateWinds &windinterp)
  {
    /* corrector velocities based on predicted coords */
    const double vel1 = windinterp.interp_uvel();

    /* predictor coords given velocity at previous coords */
    windinterp.coord1 += vel1 * delt; // move by u wind

    /* corrector velocities based on predicted coords */
    const double corrvel1 = windinterp.interp_uvel();

    /* predicted-corrected change to superdrop coords */
    const double delta1((vel1 + corrvel1) * (delt / 2));
  }

  KOKKOS_INLINE_FUNCTION
  double delta_coord2(Superdrop &drop
                      InterpolateWinds &windinterp)
  {
    /* corrector velocities based on predicted coords */
    const double vel2 = windinterp.interp_vvel();

    /* predictor coords given velocity at previous coords */
    windinterp.coord2 += vel2 * delt; // move by v wind

    /* corrector velocities based on predicted coords */
    const double corrvel2 = windinterp.interp_vvel();

    /* predicted-corrected change to superdrop coords */
    const double delta2((vel2 + corrvel2) * (delt / 2));
  }

public:
  PredCorrMotion(const unsigned int motionstep)
      : interval(motionstep) {}

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
    InterpolateWinds windinterp{gbxmaps, state, gbxindex,
                                drop.coord3, drop.coord1, drop.coord2}; // TODO
    const double delta3 = delta_coord3(drop, windinterp);
    const double delta1 = delta_coord1(drop, windinterp);
    const double delta2 = delta_coord2(drop, windinterp);

    /* CFL check on predicted change to SD coords */
    cfl_criteria(gbxmaps, gbxindex, delta3, delta1, delta2);

    /* update SD coords */
    drop.increment_coords(delta3, delta1, delta2);
  }
};

#endif // PREDCORRMOTION_HPP


