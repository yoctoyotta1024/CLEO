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

template <GridboxMaps GbxMaps>
struct PredCorrMotion
{
private:
  unsigned int interval;

  struct DeltaCoords
  /* change in coordinates calculated by predictor
  corrector method with the wind velocity obtained via 
  a simple linear interpolation. Methods follows
  equations in Grabowski et al. 2018 */
  {
    double delta3; // change in coord3
    double delta1; // change in coord1
    double delta2; // change in coord2

    DeltaCoords(const unsigned int gbxindex,
                const GbxMaps &gbxmaps,
                const State &state,
                Superdrop &drop)
    {
      // const double terminal = terminalv(drop);

      // WindsAtCoord winds{gbxmaps, state, gbxindex,
      //                    drop.coord3, drop.coord1, drop.coord2};

      // /* corrector velocities based on predicted coords */
      // const double vel3 = winds.interp_wvel() - terminal;
      // const double vel1 = winds.interp_uvel();
      // const double vel2 = winds.interp_vvel();

      // /* predictor coords given velocity at previous coords */
      // winds.coord3 += vel3 * delt; // move by w wind + terminal velocity
      // winds.coord1 += vel1 * delt; // move by u wind
      // winds.coord2 += vel2 * delt; // move by v wind

      // /* corrector velocities based on predicted coords */
      // const double corrvel3 = winds.interp_wvel() - terminal;
      // const double corrvel1 = winds.interp_uvel();
      // const double corrvel2 = winds.interp_vvel();

      // /* predicted-corrected change to superdrop coords */
      // const double delta3((vel3 + corrvel3) * (delt / 2));
      // const double delta1((vel1 + corrvel1) * (delt / 2));
      // const double delta2((vel2 + corrvel2) * (delt / 2));

      // return Deltas{delta3, delta1, delta2};
    }
  };

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
    const double delta3 = deltafucn3(gbxindex, gbxmaps, state, drop);
    const double delta1 = deltafucn3(gbxindex, gbxmaps, state, drop);
    const double delta2 = deltafucn3(gbxindex, gbxmaps, state, drop);

    /* CFL check on predicted change to SD coords */
    cfl_criteria(gbxmaps, gbxindex, delta3, delta1, delta2);

    /* update SD coords */
    drop.increment_coords(delta3, delta1, delta2);
  }
};

#endif // PREDCORRMOTION_HPP