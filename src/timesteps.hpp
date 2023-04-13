// Author: Clara Bayley
// File: timesteps.hpp
/* Structs for handling values of
timestep variables for SDM */

#ifndef TIMESTEPS_HPP
#define TIMESTEPS_HPP

#include <chrono>
#include <ratio>
#include <stdexcept>
#include <algorithm>

#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

namespace dlc = dimless_constants;

using model_step = std::chrono::duration<int, std::ratio<1, 100000>>;
using secd = std::chrono::duration<double>;

inline int realtime2timestep(const double TSTEP)
/* convert TSTEP [seconds] (a double e.g. from config file)
into a dimensionless time and then into an integer no. of
model timesteps using model_step chrono function */
{
  return std::chrono::round<model_step>(secd{TSTEP / dlc::TIME0}).count();
}

inline double timestep2realtime(const int mdlstep)
/* convert model timestep (integer) into a dimensionless time
given secd and model_step chrono functions */
{
  return secd{model_step{mdlstep}}.count() * dlc::TIME0;
}

inline double timestep2dimlesstime(const int mdlstep)
/* convert model timestep (integer) into a dimensionless time
given secd and model_step chrono functions */
{
  return secd{model_step{mdlstep}}.count();
}

struct Timesteps
{
  const int condstep;
  const int collstep;
  const int sedistep;
  const int xchangestep;
  const int outstep;
  const int tend;

  Timesteps(Config config);
  /* double's that are timesteps in config struct
  are converted into integer values of model timesteps using
  model_step and secd template functions created using std::chrono library.
  Throw error if after convertion into model timestep, any
  timestep = 0 */
};

#endif // TIMESTEPS_HPP