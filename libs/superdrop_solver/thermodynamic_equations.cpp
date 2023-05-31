// Author: Clara Bayley
// File: thermodynamic_equations.cpp
/* functions that return Left Hand
Side (LHS) of thermodynamic equations */

#include "thermodynamic_equations.hpp"

double saturation_pressure(const double temp)
/* Calculate the equilibrium vapor pressure
of water over liquid water ie. the
saturation pressure (psat). Equation taken from Bjorn Steven's
"make_tetens" python function from his module
"moist_thermodynamics.saturation_vapour_pressures"
available on gitlab. Original paper
"Murray, F. W. On the Computation of Saturation
Vapor Pressure. Journal of Applied Meteorology
and Climatology 6, 203–204 (1967)."
Note function is called with conversion
to real temp /K = T*Temp0 and from
real psat to dimensionless psat = psat/P0. */
{
  constexpr double A = 17.4146; // constants from Bjorn Gitlab originally from paper
  constexpr double B = 33.639; // ditto
  constexpr double TREF = 273.16;  // Triple point temperature [K] of water
  constexpr double PREF = 611.655; // Triple point pressure [Pa] of water

  const double T = temp * dlc::TEMP0; // real T [K]

  if (T <= 0.0)
  {
    const std::string err("psat ERROR: T must be larger than 0K. T = " +
                            std::to_string(T));
    throw std::invalid_argument(err);
  }

  return (PREF * exp(A * (T - TREF) / (T - B))) / dlc::P0; // dimensionless psat
}

double saturation_pressure_murphy_koop(const double temp)
/* Calculate the equilibrium vapor pressure
  of water over liquid water ie. the
  saturation pressure (psat). Equation taken from
  python module typhon.physics.thermodynamics.e_eq_water_mk
  with conversion to real temp /K = T*Temp0 and from
  real psat to dimensionless psat = psat/P0. */
{
  const double T = temp * dlc::TEMP0; // real T [K]

  if (T <= 0.0)
  {
    const std::string errormsg = "psat ERROR: T must be larger than 0K. T = " +
                                 std::to_string(T);
    throw std::invalid_argument(errormsg);
  }

  const double lnpsat = (54.842763 // ln(psat) [Pa]
                         - 6763.22 / T - 4.21 * log(T) + 0.000367 * T + tanh(0.0415 * (T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T));

  return exp(lnpsat) / dlc::P0; // dimensionless psat
}