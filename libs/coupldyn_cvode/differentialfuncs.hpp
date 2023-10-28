/*
 * ----- CLEO -----
 * File: differentialfuncs.hpp
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

#ifndef DIFFERENTIALFUNCS_HPP
#define DIFFERENTIALFUNCS_HPP

#include <cmath>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */

#include "cleoconstants.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/* user data structure for passing
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct pUserData
{
  size_t neq;
  realtype wmax;
  realtype tauhalf;
} * UserData;

inline double cvode_massmixingratio(const double press_vapour,
                                    const double press)
/* Calculate mass mixing ratio, qv = m_v/m_dry = rho_v/rho_dry
given the vapour pressure, pv = p_v/p_tot and total pressure p_tot */
{
  return dlc::Mr_ratio * press_vapour / (press - press_vapour);
}

inline double cvode_moistspecifcheat(const double qv, const double qc)
/* effective specific heat capacity of moist parcel
  of air (dry + water vapour + liquid water) */
{
  return dlc::Cp_dry + (dlc::Cp_v * qv) + (dlc::C_l * qc);
}

double cvode_saturationpressure(const double temp);
/* Calculate the equilibrium vapor pressure of water over liquid
water ie. the saturation pressure (psat). Equation taken from
Bjorn Steven's "make_tetens" python function from his module
"moist_thermodynamics.saturation_vapour_pressures" available
on gitlab. Original paper "Murray, F. W. On the Computation of
Saturation Vapor Pressure. Journal of Applied Meteorology
and Climatology 6, 203–204 (1967)." Note function is called
with conversion to real temp /K = T*Temp0 and from real psat
to dimensionless psat = psat/P0. */

int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data);
/* Simple function f(t,y, ydot) called by ODE solver to
  integrate ODEs over time. */


#endif // DIFFERENTIALFUNCS_HPP


