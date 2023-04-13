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
#include <string>
#include <map>

struct ThermoState
{
  const double volume;
  double time;
  
  double press;
  double temp;
  double qvap;
  double qcond;
  double wvel;
  double uvel;
  double vvel;

  ThermoState(const double vol) : volume(vol), time(), press(),
                                  temp(), qvap(), qcond(),
                                  wvel(), uvel(), vvel(){};
  inline double operator[](const int i) const;
  inline ThermoState operator-(const ThermoState &previousstate) const;
  inline bool operator==(const ThermoState &previousstate) const;
};

ThermoState ThermoState::operator-(const ThermoState &previousstate) const
{
  ThermoState delta_state(volume);

  delta_state.temp = temp - previousstate.temp;
  delta_state.qvap = qvap - previousstate.qvap;
  delta_state.qcond = qcond - previousstate.qcond;
  
  delta_state.wvel = wvel - previousstate.wvel;
  delta_state.uvel = uvel - previousstate.uvel;
  delta_state.vvel = vvel - previousstate.vvel;

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