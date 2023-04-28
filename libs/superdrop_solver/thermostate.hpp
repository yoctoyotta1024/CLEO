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
#include <utility>

struct ThermoState
{
  const double volume;
  double time;

  double press;                   // defined at centre of volume
  double temp;                    // defined at centre of volume
  double qvap;                    // defined at centre of volume
  double qcond;                   // defined at centre of volume
  std::pair<double, double> wvel; // defined on {lower, upper} z faces of volume
  std::pair<double, double> uvel; // defined on {lower, upper} x faces of volume
  std::pair<double, double> vvel; // defined on {lower, upper} y faces of volume

  ThermoState(const double vol) : volume(vol), time(), press(),
                                  temp(), qvap(), qcond(),
                                  wvel(), uvel(), vvel(){};

  ThermoState operator-(const ThermoState &prevstate) const
  {
    std::pair<double, double>
        subtract = [](const std::pair<double, double> a,
                          const std::pair<double, double> b)
    {
      return { a.first - b.first, a.second - b.second };
    };

    ThermoState delta_state(volume);
    delta_state.temp = temp - prevstate.temp;
    delta_state.qvap = qvap - prevstate.qvap;
    delta_state.qcond = qcond - prevstate.qcond;

    delta_state.wvel = subtract(wvel, prevstate.wvel);
    delta_state.uvel = subtract(uvel, prevstate.uvel);
    delta_state.vvel = subtract(vvel, prevstate.vvel);
    
    return delta_state;
  }

  bool operator==(const ThermoState &prevstate) const
  {
    if (temp == prevstate.temp &&
        qvap == prevstate.qvap &&
        qcond == prevstate.qcond)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
};

#endif // THERMOSTATE_HPP