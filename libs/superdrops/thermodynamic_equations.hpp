/*
 * ----- CLEO -----
 * File: thermodynamic_equations.hpp
 * Project: superdrops
 * Created Date: Wednesday 25th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 26th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for functions that return
 * Left Hand Side of thermodynamic equations
 * */

#ifndef THERMODYNAMIC_EQUATIONS_HPP
#define THERMODYNAMIC_EQUATIONS_HPP

#include <cmath>
// #include <stdexcept>
#include <cassert>
#include <string>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;

KOKKOS_INLINE_FUNCTION
double moist_specifc_heat(const double qvap, const double qcond)
/* (dimensionless) specific heat capacity of a moist parcel of air */
{
  return dlc::Cp_dry + dlc::Cp_v * (qvap) + dlc::C_l * (qcond);
}

KOKKOS_INLINE_FUNCTION
double supersaturation_ratio(const double press,
                             const double qvap,
                             const double psat)
/* calculate the superaturation ratio given the saturaion
pressure, psat, the ambient pressure, press, and the vapour
mass mixing ratio, qvap. supersaturation ratio = 's_ratio', 
s_ratio = p_vapour/psat (ie. relative humidity) */
{
  return (press * qvap) / ((dlc::Mr_ratio + qvap) * psat);
}

KOKKOS_FUNCTION
double saturation_pressure(const double temp);
/* Calculate the equilibrium vapor pressure of water over liquid
water ie. the saturation pressure (psat). Equation taken from
Bjorn Steven's "make_tetens" python function from his module
"moist_thermodynamics.saturation_vapour_pressures" available
on gitlab. Original paper "Murray, F. W. On the Computation of
Saturation Vapor Pressure. Journal of Applied Meteorology
and Climatology 6, 203–204 (1967)." Note function is called
with conversion to real temp /K = T*Temp0 and from real psat
to dimensionless psat = psat/P0. */

KOKKOS_FUNCTION
double saturation_pressure_murphy_koop(const double temp);
/* Calculate the equilibrium vapor pressure of water over
liquid water ie. the saturation pressure (psat). Equation taken from
python module typhon.physics.thermodynamics.e_eq_water_mk
with conversion to real temp /K = T*Temp0 and from
real psat to dimensionless psat = psat/P0. */



/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION
double saturation_pressure(const double temp)
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

KOKKOS_FUNCTION
double saturation_pressure_murphy_koop(const double temp)
/* Calculate the equilibrium vapor pressure of water over
liquid water ie. the saturation pressure (psat). Equation taken
from python module typhon.physics.thermodynamics.e_eq_water_mk
with conversion to real temp /K = T*Temp0 and from
real psat to dimensionless psat = psat/P0. */
{
  assert((temp > 0) && "psat ERROR: temperature must be larger than 0K.");

  const double T(temp * dlc::TEMP0); // real T [K]

  const double lnpsat = (54.842763 // ln(psat) [Pa]
                         - 6763.22 / T - 4.21 * log(T) +
                         0.000367 * T +
                         tanh(0.0415 * (T - 218.8)) *
                             (53.878 - 1331.22 / T -
                              9.44523 * log(T) + 0.014025 * T));

  return std::exp(lnpsat) / dlc::P0; // dimensionless psat
}

#endif // THERMODYNAMIC_EQUATIONS_HPP
