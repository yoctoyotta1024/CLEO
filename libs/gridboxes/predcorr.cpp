/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: predcorr.cpp
 * Project: gridboxes
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for change in a superdroplet's coords
 * using predictor-corrector method for motion of droplet
 * given a formula for its terminal velocity and the wind
 * velocity obtained via a simple linear interpolation.
 * Methods follows equations in Grabowski et al. 2018
 */

#include "./predcorr.hpp"

/* Given [X = z,x or y] wind velocity component, vel, that is
defined on the faces of a gridbox at {lower, upper} [X] bounds,
return wind at [X] coord. Method is 'simple' linear interpolation
from Grabowski et al. (2018). coord use in interpolation is
limited to lower_bound <= coord <= upper_bound. */
KOKKOS_FUNCTION
double interpolation(const Kokkos::pair<double, double> bounds,
                     const Kokkos::pair<double, double> vel, const double sdcoord) {
  const auto coord = double{Kokkos::fmin(
      bounds.second, Kokkos::fmax(bounds.first, sdcoord))};  // limit coord to within bounds

  const auto alpha = double{(coord - bounds.first) / (bounds.second - bounds.first)};

  const auto interp =
      double{alpha * vel.second + (1 - alpha) * vel.first};  // simple linear interpolation

  return interp;
}
