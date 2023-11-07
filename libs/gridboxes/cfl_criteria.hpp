/*
 * ----- CLEO -----
 * File: cfl_criteria.hpp
 * Project: gridboxes
 * Created Date: Tuesday 7th November 2023
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
 */

#ifndef CFL_CRITERIA_HPP
#define CFL_CRITERIA_HPP

#include "./gridboxmaps.hpp"

inline bool cfl_criterion(const double gridstep,
                          const double sdstep)
/* sdstep = change in superdroplet coordinate position.
returns *false* if cfl criterion, C = sdstep / gridstep, > 1 */
{
  return (std::abs(sdstep) <= std::abs(gridstep));
}

template <GridboxMaps GbxMaps>
bool cfl_criteria(const Maps4GridBoxes &gbxmaps,
                  const unsigned int gbxindex,
                  const double delta3,
                  const double delta1,
                  const double delta2)
/* returns false if any of z, x or y (3,1,2) directions
  do not meet their cfl criterion. For each direction,
  Criterion is C = delta[X] / gridstep =< 1 where the
  gridstep is calculated from the gridbox boundaries map */
{
  double gridstep(gbxmaps.coord3bounds(gbxindex).second -
                  gbxmaps.coord3bounds(gbxindex).first);
  bool cfl(cfl_criterion(gridstep, delta3));

  gridstep = gbxmaps.coord1bounds(gbxindex).second -
             gbxmaps.coord1bounds(gbxindex).first;
  cfl = (cfl_criterion(gridstep, delta1) && cfl);

  gridstep = gbxmaps.coord2bounds(gbxindex).second -
             gbxmaps.coord2bounds(gbxindex).first;
  cfl = (cfl_criterion(gridstep, delta2) && cfl);

  if (!cfl)
  {  
    throw std::invalid_argument("CFL criteria for SD motion not met."
                                "Consider reducing sdmotion timestep");
  }

  return cfl;
}

#endif // CFL_CRITERIA_HPP