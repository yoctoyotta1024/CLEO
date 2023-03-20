// Author: Clara Bayley
// File: thermostate.hpp
/* Header for objects used in
handling of thermodynamic variables
(pressure, temperaure, q_vapour,
q_condensate, volume, time) of SDM */

#ifndef THERMOSTATE_HPP
#define THERMOSTATE_HPP

#include <iostream>
#include <stdexcept>

struct ThermoState
{
  double press;
  double temp;
  double qvap;
  double qcond;
  double volume;
  double time;

  inline ThermoState();
  inline double operator[](const int i) const;
  inline ThermoState operator-(const ThermoState &previousstate) const;
  inline bool operator==(const ThermoState &previousstate) const;
};

ThermoState::ThermoState()
    : press(), temp(), qvap(),
      qcond(), volume(), time()
{
}

double ThermoState::operator[](const int i) const
{
  if (i == -1)
  {
    return time;
  }
  else if (i == 0)
  {
    return press;
  }
  else if (i == 1)
  {
    return temp;
  }
  else if (i == 2)
  {
    return qvap;
  }
  else if (i == 3)
  {
    return qcond;
  }
  else if (i == 4)
  {
    return volume;
  }
  else
  {
    const std::string errormsg = "index out of range for thermo state";
    throw std::invalid_argument(errormsg);
  }
}

ThermoState ThermoState::operator-(const ThermoState &previousstate) const
{
  ThermoState delta_state;

  delta_state.temp = temp - previousstate.temp;
  delta_state.qvap = qvap - previousstate.qvap;
  delta_state.qcond = qcond - previousstate.qcond;

  return delta_state;
}

bool ThermoState::operator==(const ThermoState &previousstate) const
{
  if (temp == previousstate.temp &&
      qvap == previousstate.qvap &&
      qcond == previousstate.qcond)
  {
    return true;
  }
  else
  {
    return false;
  }
}
#endif // THERMOSTATE_HPP