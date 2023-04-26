// Author: Clara Bayley
// File: "sdmotion.cpp"
/* Implementation of some functions
related to moving superdroplets by
updating their coordinates according
to equations of motion */

#include "sdmotion.hpp"

bool cfl_criteria(const Maps4GridBoxes &gbxmaps,
                  const unsigned int gbxindex,
                  const double delta3, const double delta1,
                  const double delta2)
/* returns false if any of z, x or y (3,1,2) directions
  do not meet their cfl criterion. For each direction,
  Criterion is C = delta[X] / gridstep =< 1 where the
  gridstep is calculated from the gridbox boundaries map */
{
  double gridstep(gbxmaps.get_bounds_z(gbxindex).second -
                  gbxmaps.get_bounds_z(gbxindex).first);
  bool cfl(cfl_criterion(gridstep, delta3));

  gridstep = gbxmaps.get_bounds_x(gbxindex).second -
             gbxmaps.get_bounds_x(gbxindex).first;
  cfl = cfl_criterion(gridstep, delta1);

  gridstep = gbxmaps.get_bounds_y(gbxindex).second -
             gbxmaps.get_bounds_y(gbxindex).first;
  cfl = cfl_criterion(gridstep, delta2);

  if (!cfl)
  {  
    throw std::invalid_argument("CFL criteria for SD motion not met."
                                "Consider reducing sdmotion timestep");
  }

  return cfl;
}

void NoInterpMoveWithSedimentation::
    change_superdroplet_coords(const Maps4GridBoxes &gbxmaps,
                               const GridBox &gbx,
                               Superdrop &drop) const
/* very crude method to forward timestep the velocity
using the velocity from the gridbox thermostate, ie.
without interpolation to the SD position and using
single step forward euler method to integrate dx/dt */
{
  const double delta3 = deltacoord(gbx.state.wvel - terminalv(drop)); // w wind + terminal velocity
  const double delta1 = deltacoord(gbx.state.uvel); // u component of wind velocity
  const double delta2 = deltacoord(gbx.state.vvel); // v component of wind velocity (y=2)

  cfl_criteria(gbxmaps, gbx.gbxindex, delta3, delta1, delta2);

  drop.coord3 += delta3;
  drop.coord1 += delta1;
  drop.coord2 += delta2;
}

void MoveWith2DFixedFlow::
    change_superdroplet_coords(const Maps4GridBoxes &gbxmaps,
                               const GridBox &gbx,
                               Superdrop &drop) const
/* Use predictor-corrector scheme from Grabowksi et al. 2018
(similar to Arabas et al. 2015) to update a SD position.
The velocity required for this scheme is determined
from the PrescribedFlow2D instance */
{
  auto deltas = predictor_corrector(drop.coord3, drop.coord1);
  
  cfl_criteria(gbxmaps, gbx.gbxindex, deltas.first, deltas.second, 0.0);

  drop.coord3 += deltas.first;
  drop.coord1 += deltas.second; 
}

std::pair<double, double>
MoveWith2DFixedFlow::predictor_corrector(const double coord3,
                                         const double coord1) const
/* returns change in (z,x) coordinates = (delta3, delta1)
obtained using predictor-corrector method and velocities
calculated from a Prescribed2DFlow */
{ 
  const double vel3 = flow2d.prescribed_wvel(coord3, coord1); // w wind from prescribed 2D flow
  const double vel1 = flow2d.prescribed_uvel(coord3, coord1); // u wind from prescribed 2D flow

  const double pred3(coord3 + vel3 * delt);
  const double pred1(coord1 + vel1 * delt);

  const double corrvel3 = flow2d.prescribed_wvel(pred3, pred1); // w wind from prescribed 2D flow
  const double corrvel1 = flow2d.prescribed_uvel(pred3, pred1); // w wind from prescribed 2D flow

  const double delta3((vel3 + corrvel3) * (delt / 2));
  const double delta1((vel1 + corrvel1) * (delt / 2));

  return std::pair<double, double>(delta3, delta1);
}