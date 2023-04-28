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

double WindsAtCoord::interpolate_wind(const std::pair<double, double> bounds,
                                      const std::pair<double, double> vel,
                                      const double coord) const
/* Given [X = z,x or y] wind velocity component, vel, that is
defined on the faces of a gridbox at {lower, upper} [X] bounds,
return wind at [X] coord. Method is 'simple' linear interpolation
from Grabowski et al. (2018) */
{
  const double alpha((coord - bounds.first) / (bounds.second - bounds.first));

  const double interpolated_vel(alpha*vel.second + (1-alpha)*vel.first);

  return interpolated_vel;
} 

double WindsAtCoord::interp_wvel() const
/* returns w wind velocity at z=coord3 for gridbox gbxindex */
{
  return interpolate_wind(gbxmaps.get_bounds_z(gbxindex),
                          state.wvel, coord3);
}

double WindsAtCoord::interp_uvel() const
/* returns u wind velocity at x=coord1 for gridbox gbxindex */
{
  return interpolate_wind(gbxmaps.get_bounds_x(gbxindex),
                          state.uvel, coord1); 
}

double WindsAtCoord::interp_vvel() const
/* returns v wind velocity at y=coord2 for gridbox gbxindex */
{
  return interpolate_wind(gbxmaps.get_bounds_y(gbxindex),
                          state.vvel, coord2); 
}