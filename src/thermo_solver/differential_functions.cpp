// Author: Clara Bayley
// File: differential_functions.cpp
/* ODE functions which are solved by
 CVODE ode solver to model evolution of the
thermodyanmics (p, temp, qv and qc) over time */

#include "differential_functions.hpp"

static double dp_dt(const double t, const double wmax, const double tauhalf);
static double moist_specifc_heat(const double qv, const double qc);
static double dtemp_dt_adia(const int k, const double pdot, const N_Vector &y);

int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data)
/* Simple function f(t,y, ydot) called by ODE solver to
  integrate ODEs over time. */
{

  constexpr int NVARS = 4; // no. of (distinct) variables per grid box

  UserData data = (UserData)user_data;
  const size_t neq = data->neq;
  const bool doThermo = data->doThermo;
  const double wmax = data->wmax;
  const double tauhalf = data->tauhalf;

  for (size_t k = 0; k < neq; k += NVARS) // loop over grid boxes
  {
    if (doThermo)
    {
      double pdot;
      pdot = dp_dt(t, wmax, tauhalf);
      NV_Ith_S(ydot, k) = pdot;
      NV_Ith_S(ydot, k + 1) = dtemp_dt_adia(k, pdot, y);
    }
    else
    {
      NV_Ith_S(ydot, k) = 0;
      NV_Ith_S(ydot, k + 1) = 0;
    }

    NV_Ith_S(ydot, k + 2) = 0;
    NV_Ith_S(ydot, k + 3) = 0;
  }

  return 0;
}

static double dp_dt(const double t, const double wmax, const double tauhalf)
/* dp/dt differential equation (dimensionless)
  describing pressure evolution over time.
  note: true dP/dt = dp/dt * P0/TIME0 */
{
  double pdot, profile, w, z;

  constexpr double zg = 0.0 / (dlc::W0 * dlc::TIME0);                    // dimensionless z value at ground level
  constexpr double tempg = 273.15 / dlc::TEMP0;                          // dimensionless temperature at zg
  constexpr double pg = 100000.0 / dlc::P0;                              // dimensionless pressure at zg
  constexpr double lpsrate = 0.0062 / dlc::TEMP0 * dlc::W0 * dlc::TIME0; // dimensionless moist adiabatic lapse rate
  constexpr double gamma = DC::G / (DC::RGAS_DRY * 0.0062) - 1.0;        // constant in dry adiabatic expansion
  constexpr double dp_dt_const = -dlc::W0 * dlc::TIME0 * DC::G / (DC::RGAS_DRY * dlc::TEMP0) * pg / tempg;

  /* comment this in for sinusoidally time dependent velocity, w */
  w = wmax * sin(t / tauhalf);                 // sinusoidal velocity profile
  z = wmax * tauhalf * (1 - cos(t / tauhalf)); // sinusoidal z coordinate

  /* or comment this in for constant velocity, w*/
  // w = wmax * 2/ M_PI;                                 // non-sinusoidal, constant velocity profile
  // z = w * t;                                          // linear z coordinate

  profile = 1.0 - lpsrate / tempg * (z - zg); // characteristic function for pressure profile as
  profile = pow(profile, gamma);              //      a funciton of time (ie. height via z=w*t)

  pdot = dp_dt_const * profile * w;

  return pdot;
}

static double moist_specifc_heat(const double qv, const double qc)
/* effective specific heat capacity of moist parcel
  of air (dry + water vapour + liquid water) */
{
  return dlc::Cp_dry + dlc::Cp_v * qv + dlc::C_l * qc;
}

static double dtemp_dt_adia(const int k, const double pdot, const N_Vector &y)
/* dtemp/dt differential equation describing
  temperature evolution solely due to pressure
  changes in parcel for adiabatic process (no heat loss).
  Parcel has water vapour mass mixing ratio (m_v/m_dry) = qv and
  liquid water mass mixing ratio (m_c/m_dry) = qc.
  note: True dTemp/dt = dtemp * TEMP0/TIME0  */
{
  double tempdot, rho_d, cp_m;
  double p, temp, qv, qc;

  p = NV_Ith_S(y, k);
  temp = NV_Ith_S(y, k + 1);
  qv = NV_Ith_S(y, k + 2);
  qc = NV_Ith_S(y, k + 3);

  rho_d = dlc::Mr_ratio / (dlc::Mr_ratio + qv) * p / temp; // density of dry parcel (p_dry/temp)

  cp_m = moist_specifc_heat(qv, qc); // moist specific heat capacity

  tempdot = dlc::Rgas_dry / (rho_d * cp_m) * pdot;

  return tempdot;
}
