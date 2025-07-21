/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: differentialfuncs.cpp
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
 * thermodynamics (p, temp, qv and qc) over time
 */

#include "coupldyn_cvode/differentialfuncs.hpp"

/* Calculate the equilibrium vapor pressure of water over liquid
water ie. the saturation pressure (psat). Equation taken from
Bjorn Steven's "make_tetens" python function from his module
"moist_thermodynamics.saturation_vapour_pressures" available
on gitlab. Original paper "Murray, F. W. On the Computation of
Saturation Vapor Pressure. Journal of Applied Meteorology
and Climatology 6, 203–204 (1967)." Note function is called
with conversion to real temp /K = T*Temp0 and from real psat
to dimensionless psat = psat/P0. */
double cvode_saturationpressure(const double temp) {
  assert((temp > 0) && "psat ERROR: temperature must be larger than 0K.");

  constexpr double A = 17.4146;     // constants from Bjorn Gitlab originally from paper
  constexpr double B = 33.639;      // ditto
  constexpr double TREF = 273.16;   // Triple point temperature [K] of water
  constexpr double PREF = 611.655;  // Triple point pressure [Pa] of water

  const auto T = double{temp * dlc::TEMP0};  // real T [K]

  return (PREF * std::exp(A * (T - TREF) / (T - B))) / dlc::P0;  // dimensionless psat
}

/* dp/dt differential equation (dimensionless)
  describing pressure evolution over time.
  _Note:_ true dP/dt = dp/dt * P0/TIME0 */
static double dp_dt(const double t, const double wmax, const double tauhalf) {
  constexpr double zg = 0.0 / (dlc::W0 * dlc::TIME0);  // dimensionless z value at ground level
  constexpr double tempg = 273.15 / dlc::TEMP0;        // dimensionless temperature at zg
  constexpr double pg = 100000.0 / dlc::P0;            // dimensionless pressure at zg
  constexpr double lpsrate =
      0.0062 / dlc::TEMP0 * dlc::W0 * dlc::TIME0;  // dimensionless moist adiabatic lapse rate
  constexpr double gammafac =
      DC::G / (DC::RGAS_DRY * 0.0062) - 1.0;  // constant in dry adiabatic expansion
  constexpr double dp_dt_const =
      -dlc::W0 * dlc::TIME0 * DC::G / (DC::RGAS_DRY * dlc::TEMP0) * pg / tempg;

  /* comment this in for sinusoidally time dependent velocity, w */
  double w(wmax * sin(t / tauhalf));                  // sinusoidal velocity profile
  double z(wmax * tauhalf * (1 - cos(t / tauhalf)));  // sinusoidal z coordinate

  /* or comment this in for constant velocity, w*/
  // double w(wmax * 2/ M_PI);                                 // non-sinusoidal, constant velocity
  // profile double z(w * t);                                          // linear z coordinate

  double profile(1.0 -
                 lpsrate / tempg * (z - zg));  // characteristic function for pressure profile as
  profile = std::pow(profile, gammafac);       //      a funciton of time (ie. height via z=w*t)

  double pdot(dp_dt_const * profile * w);
  return pdot;
}

/* dtemp/dt differential equation describing
  temperature evolution solely due to pressure
  changes in parcel for adiabatic process (no heat loss).
  Parcel has water vapour mass mixing ratio (m_v/m_dry) = qv and
  liquid water mass mixing ratio (m_c/m_dry) = qc.
  _Note:_ True dTemp/dt = dtemp * TEMP0/TIME0  */
static double dtemp_dt_adia(const int k, const double pdot, const N_Vector &y) {
  double p(NV_Ith_S(y, k));
  double temp(NV_Ith_S(y, k + 1));
  double qv(NV_Ith_S(y, k + 2));
  double qc(NV_Ith_S(y, k + 3));

  double rho_d(dlc::Mr_ratio / (dlc::Mr_ratio + qv) * p /
               temp);  // density of dry parcel (p_dry/temp)

  double cp_m(cvode_moistspecifcheat(qv, qc));  // moist specific heat capacity

  double tempdot(dlc::Rgas_dry / (rho_d * cp_m) * pdot);

  return tempdot;
}

/* Simple function f(t,y, ydot) called by ODE solver to
  integrate ODEs over time. */
int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  constexpr int NVARS = 4;  // no. of (distinct) variables per grid box

  UserData data = (UserData)user_data;
  const size_t neq(data->neq);
  const auto wmax = double{data->wmax};
  const auto tauhalf = double{data->tauhalf};

  // loop over grid boxes
  for (size_t k = 0; k < neq; k += NVARS) {
    double pdot(dp_dt(t, wmax, tauhalf));
    NV_Ith_S(ydot, k) = pdot;
    NV_Ith_S(ydot, k + 1) = dtemp_dt_adia(k, pdot, y);
    NV_Ith_S(ydot, k + 2) = 0;
    NV_Ith_S(ydot, k + 3) = 0;
  }

  return 0;
}
