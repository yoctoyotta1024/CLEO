/*
 * ----- CLEO -----
 * File: predcorrmotion.cpp
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
 * Functionality for motion using predictor-corrector
 * method to update a superdroplet's coordinates given
 * a formula for its terminal velocity and the wind
 * velocity obtained via a simple linear interpolation.
 * Methods follows equations in Grabowski et al. 2018
 */

#include "./predcorrmotion.hpp"

/* method to interpolate (w, u, and v) wind velocity components
defined on (coord3, coord1 and coord2) faces of a gridbox
to a superdroplet's coordinates at (coord3, coord1, coord2) */

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