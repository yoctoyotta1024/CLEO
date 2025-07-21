/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: differentialfuncs.hpp
 * Project: coupldyn_cvode
 * Created Date: Saturday 28th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for ODEs which are solved by
 * CVODE ode solver to model evolution of the
 * thermodynamics (p, temp, qv and qc) over time */

#ifndef LIBS_COUPLDYN_CVODE_DIFFERENTIALFUNCS_HPP_
#define LIBS_COUPLDYN_CVODE_DIFFERENTIALFUNCS_HPP_

#include <nvector/nvector_serial.h> /* access to serial N_Vector            */

#include <cassert>
#include <cmath>

#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/* user data structure for passing
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct pUserData {
  size_t neq;
  realtype wmax;
  realtype tauhalf;
} *UserData;

/* Calculate mass mixing ratio, qv = m_v/m_dry = rho_v/rho_dry
given the vapour pressure, pv = p_v/p_tot and total pressure p_tot */
inline double cvode_massmixingratio(const double press_vapour, const double press) {
  return dlc::Mr_ratio * press_vapour / (press - press_vapour);
}

/* effective specific heat capacity of moist parcel
  of air (dry + water vapour + liquid water) */
inline double cvode_moistspecifcheat(const double qv, const double qc) {
  return dlc::Cp_dry + (dlc::Cp_v * qv) + (dlc::C_l * qc);
}

/* Calculate the equilibrium vapor pressure of water over liquid
water ie. the saturation pressure (psat). Equation taken from
Bjorn Steven's "make_tetens" python function from his module
"moist_thermodynamics.saturation_vapour_pressures" available
on gitlab. Original paper "Murray, F. W. On the Computation of
Saturation Vapor Pressure. Journal of Applied Meteorology
and Climatology 6, 203–204 (1967)." Note function is called
with conversion to real temp /K = T*Temp0 and from real psat
to dimensionless psat = psat/P0. */
double cvode_saturationpressure(const double temp);

/* Simple function f(t,y, ydot) called by ODE solver to
  integrate ODEs over time. */
int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#endif  // LIBS_COUPLDYN_CVODE_DIFFERENTIALFUNCS_HPP_
