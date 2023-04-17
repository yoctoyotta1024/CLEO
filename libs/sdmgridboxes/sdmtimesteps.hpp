// Author: Clara Bayley
// File: sdmtimesteps.hpp
/* Structs for handling values of
timestep variables for SDM */

#ifndef TIMESTEPS_HPP
#define TIMESTEPS_HPP

#include <chrono>
#include <ratio>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "../claras_SDconstants.hpp"

namespace dlc = dimless_constants;

using model_step = std::chrono::duration<int, std::ratio<1, 100000>>;
using secd = std::chrono::duration<double>;

inline int realtime2step(const double TSTEP)
/* convert TSTEP [seconds] (e.g. a double from Config struct)
into a dimensionless time and then into an integer no. of
model steps using model_step chrono function */
{
  return std::chrono::round<model_step>(secd{TSTEP / dlc::TIME0}).count();
}

inline double step2realtime(const int mdlstep)
/* convert model step (integer) into a time [seconds]
given secd and model_step chrono functions */
{
  return secd{model_step{mdlstep}}.count() * dlc::TIME0;
}

inline double step2dimlesstime(const int mdlstep)
/* convert model timestep (integer) into a dimensionless
time given secd and model_step chrono functions */
{
  return secd{model_step{mdlstep}}.count();
}

struct SDMTimesteps
/* integer intervals (timesteps) involved
in running coupled model */
{
  const int condsubstep;
  const int collsubstep;
  const int motionstep;
  const int couplstep;
  const int t_end;

  SDMTimesteps(const double CONDTSTEP, const double COLLTSTEP,
            const double MOTIONTSTEP, const double COUPLTSTEP,
            const double T_END);
  /* (dimensionless) double's that are timesteps in config struct
  are converted into integer values of model timesteps using
  model_step and secd template functions created using std::chrono library.
  Throw error if after convertion into model timestep, any
  timestep = 0 */
};

#endif // TIMESTEPS_HPP