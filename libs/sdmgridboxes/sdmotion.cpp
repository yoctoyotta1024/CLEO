// Author: Clara Bayley
// File: "sdmotion.cpp"
/* Implementation of some functions
related to moving superdroplets by
updating their coordinates according
to equations of motion */

#include "sdmotion.hpp"

bool cfl_criteria(const Maps4GridBoxes &gbxmaps,
                  const unsigned int gbxindex,
                  const double delt, const double wvel,
                  const double uvel, const double vvel)
/* returns false if any of z,x or y directions
  do not meet their cfl criterion. For each direction,
  Criterion is C = velocity_component*delt / gridstep =< 1
  where the gridstep is calculated from the gridbox boundaries
  map (in the same direction as the velocity component) */
{
  double gridstep(gbxmaps.get_bounds_z(gbxindex).second -
                  gbxmaps.get_bounds_z(gbxindex).first);
  bool cfl(cfl_criterion(gridstep, wvel, delt));

  gridstep = gbxmaps.get_bounds_x(gbxindex).second -
             gbxmaps.get_bounds_x(gbxindex).first;
  cfl = cfl_criterion(gridstep, uvel, delt);

  gridstep = gbxmaps.get_bounds_y(gbxindex).second -
             gbxmaps.get_bounds_y(gbxindex).first;
  cfl_criterion(gridstep, vvel, delt);

  return cfl;
}

void NoInterpMoveWithSedimentation::
    change_superdroplet_coords(const Maps4GridBoxes &gbxmaps,
                               const GridBox &gbx,
                               Superdrop &drop) const
/* very crude method to forward timestep the velocity
using the velocity from the gridbox thermostate, ie.
without interpolation to the SD position and using
forward euler method to integrate dx/dt instead of
a better method e.g. a predictor-corrector scheme */
{
  const double vel3 = gbx.state.wvel - terminalv(drop); // w wind + terminal velocity
  const double vel1 = gbx.state.uvel;                   // u component of wind velocity
  const double vel2 = gbx.state.vvel;                   // v component of wind velocity (y=2)

  const bool cfl = cfl_criteria(gbxmaps, gbx.gbxindex, delt,
                                vel3, vel1, vel2);

  if (cfl)
  {
    drop.coord3 += deltacoord(vel3);
    drop.coord1 += deltacoord(vel1);
    drop.coord2 += deltacoord(vel2);
  }
  else
  {
    throw std::invalid_argument("CFL criteria for SD motion not met."
                                "Consider reducing sdmotion timestep");
  }
}

void MoveWith2DFixedFlow::
    change_superdroplet_coords(const Maps4GridBoxes &gbxmaps,
                               const GridBox &gbx,
                               Superdrop &drop) const
/* Use predictor-corrector scheme from Grabowksi et al. 2018
(similar to Arabas et al. 2015) to update a SD position.
The velocity required for this scheme is determined
from the fixed 2D flow wiht constant density from
Arabas et al. 2015 with lengthscales X = 2*pi*xscale,
and Z = pi*zscale. No change to y component (vel2=0.0). */
{
  const double vel3 = flow2d.prescribed_wvel(drop.coord3, drop.coord1); // w wind from prescribed 2D flow
  const double vel1 = flow2d.prescribed_uvel(drop.coord3, drop.coord1); // u wind from prescribed 2D flow

  const bool cfl = cfl_criteria(gbxmaps, gbx.gbxindex, delt,
                                vel3, vel1, 0.0);

  if (cfl)
  {
    drop.coord3 += deltacoord(vel3);
    drop.coord1 += deltacoord(vel1);
  }
  else
  {
    throw std::invalid_argument("CFL criteria for SD motion not met."
                                "Consider reducing sdmotion timestep");
  }
}