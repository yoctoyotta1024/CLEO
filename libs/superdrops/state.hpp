/*
 * ----- CLEO -----
 * File: state.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for functions and structures related to detectors
 * which track data for output (e.g. of microphysical processes)
 * in gridboxes
 */


#ifndef STATE_HPP 
#define STATE_HPP 

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

struct State
/* variables which define the state of a certain volume
(e.g. inside a gridbox) at a given time. State variables
are for example thermodynamic variables required for
CLEO SDM such as temperature, pressure etc. */
{
private:
  double volume;

public:
  double press;                   // defined at centre of volume
  double temp;                    // defined at centre of volume
  double qvap;                    // defined at centre of volume
  double qcond;                   // defined at centre of volume
  Kokkos::pair<double, double> wvel; // defined on {lower, upper} z faces of volume
  Kokkos::pair<double, double> uvel; // defined on {lower, upper} x faces of volume
  Kokkos::pair<double, double> vvel; // defined on {lower, upper} y faces of volume

  KOKKOS_INLINE_FUNCTION State() = default; // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~State() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION State(const double volume,
                               const double press, const double temp,
                               const double qvap, const double qcond,
                               const Kokkos::pair<double, double> wvel,
                               const Kokkos::pair<double, double> uvel,
                               const Kokkos::pair<double, double> vvel)
      : volume(volume),
        press(press), temp(temp),
        qvap(qvap), qcond(qcond),
        wvel(wvel), uvel(uvel), vvel(vvel) {}

  KOKKOS_INLINE_FUNCTION
  auto get_volume() const { return volume; }

  KOKKOS_INLINE_FUNCTION
  double wvelcentre() const
  /* return wvel defined at centre of volume */
  {
    return (wvel.first + wvel.second) / 2;
  }

  KOKKOS_INLINE_FUNCTION
  double uvelcentre() const
  /* return uvel defined at centre of volume */
  {
    return (uvel.first + uvel.second) / 2;
  }

  KOKKOS_INLINE_FUNCTION
  double vvelcentre() const
  /* return vvel defined at centre of volume */
  {
    return (vvel.first + vvel.second) / 2;
  }
};

#endif // STATE_HPP 