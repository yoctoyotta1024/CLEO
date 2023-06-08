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
  cfl = (cfl_criterion(gridstep, delta1) && cfl);

  gridstep = gbxmaps.get_bounds_y(gbxindex).second -
             gbxmaps.get_bounds_y(gbxindex).first;
  cfl = (cfl_criterion(gridstep, delta2) && cfl);

  if (!cfl)
  {  
    throw std::invalid_argument("CFL criteria for SD motion not met."
                                "Consider reducing sdmotion timestep");
  }

  return cfl;
}