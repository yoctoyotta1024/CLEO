// Author: Clara Bayley
// File: thermodynamic_equations.hpp
/* Header file for functions that return
Left Hand Side of thermodynamic equations */

#ifndef THERMODYNAMIC_EQUATIONS_HPP
#define THERMODYNAMIC_EQUATIONS_HPP

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>

#include "../claras_SDconstants.hpp"

namespace dlc = dimless_constants;

inline double vapourpressure_2_massmixratio(const double press_vapour,
                                            const double press)
/* Calculate mass mixing ratio, qv = m_v/m_dry = rho_v/rho_dry
given the vapour pressure, pv = p_v/p_tot. */
{
  return dlc::Mr_ratio * press_vapour / (press - press_vapour);
}

inline double moist_specifc_heat(const double qvap, const double qcond)
/* (dimensionless) specific heat capacity of a moist parcel of air */
{
  return dlc::Cp_dry + dlc::Cp_v * (qvap) + dlc::C_l * (qcond);
}

inline double supersaturation_ratio(const double press, const double qvap,
                                    const double psat)
/* calculate the superaturation ratio given the saturaion pressure, psat,
the ambient pressure, press, and the vapour mass mixing ratio, qvap.
supersaturation ratio = s_ratio = p_vapour/psat (ie. relative humidity) */
{
  return (press * qvap) / ((dlc::Mr_ratio + qvap) * psat);
}

double saturation_pressure(const double temp);
/* Calculate the equilibrium vapor pressure of water over
liquid water ie. the saturation pressure (psat).
Equation taken from Bjorn Steven's "make_tetens" python function
from his module "moist_thermodynamics.saturation_vapour_pressures".
function called with conversion to real temp /K = T*Temp0 and from
real psat to dimensionless psat = psat/P0. */

double saturation_pressure_murphy_koop(const double temp);
/* Calculate the equilibrium vapor pressure of water over
liquid water ie. the saturation pressure (psat). Equation taken from
python module typhon.physics.thermodynamics.e_eq_water_mk
with conversion to real temp /K = T*Temp0 and from
real psat to dimensionless psat = psat/P0. */

#endif // THERMODYNAMIC_EQUATIONS_HPP