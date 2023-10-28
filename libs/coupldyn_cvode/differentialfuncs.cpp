/*
 * ----- CLEO -----
 * File: differentialfuncs.cpp
 * Project: coupldyn_cvode
 * Created Date: Saturday 28th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 28th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for ODEs which are solved by
 * CVODE ode solver to model evolution of the
 * thermodynamics (p, temp, qv and qc) over time */


#include "./differentialfuncs.hpp"

double cvode_saturationpressure(const double temp)
/* Calculate the equilibrium vapor pressure of water over liquid
water ie. the saturation pressure (psat). Equation taken from
Bjorn Steven's "make_tetens" python function from his module
"moist_thermodynamics.saturation_vapour_pressures" available
on gitlab. Original paper "Murray, F. W. On the Computation of
Saturation Vapor Pressure. Journal of Applied Meteorology
and Climatology 6, 203–204 (1967)." Note function is called
with conversion to real temp /K = T*Temp0 and from real psat
to dimensionless psat = psat/P0. */
{
  assert((temp > 0) && "psat ERROR: temperature must be larger than 0K.");

  constexpr double A = 17.4146; // constants from Bjorn Gitlab originally from paper
  constexpr double B = 33.639; // ditto
  constexpr double TREF = 273.16;  // Triple point temperature [K] of water
  constexpr double PREF = 611.655; // Triple point pressure [Pa] of water

  const double T(temp * dlc::TEMP0); // real T [K]

  return (PREF * std::exp(A * (T - TREF) / (T - B))) / dlc::P0; // dimensionless psat
}