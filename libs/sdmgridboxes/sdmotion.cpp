// Author: Clara Bayley
// File: "sdmotion.cpp"
/* Implementation of some functions
related to moving superdroplets by
updating their coordinates according
to equations of motion */

#include "sdmotion.hpp"

bool cfl_criteria(const Maps4GridBoxes &gbxmaps,
                  const unsigned int gbxindex,
                  const double delt,const double wvel,
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