// Author: Clara Bayley
// File: timesteps.cpp
/* Structs for handling values of
timstep variables for SDM */

#include "timesteps.hpp"

ModelTimesteps::ModelTimesteps(const double CONDTSTEP, const double COLLTSTEP,
const double MOTIONTSTEP, const double COUPLTSTEP, const double T_END)
/* (dimensionless) double's that are timesteps in config struct
are converted into integer values of model timesteps using
model_step and secd template functions created using std::chrono library.
Throw error if after convertion into model timestep, any
timestep = 0 */
    : condsubstep(realtime2step(CONDTSTEP)),
      collsubstep(realtime2step(COLLTSTEP)),
      motionstep(realtime2step(MOTIONTSTEP)),
      couplstep(realtime2step(COUPLTSTEP)),
      t_end(realtime2step(T_END))
{
  if ((condsubstep == 0) | (collsubstep == 0) | (motionstep == 0) |
      (couplstep == 0) | (t_end == 0))
  {
    const std::string err("A model step = 0, possibly due to bad conversion"
                          " of a real timestep [s]. Consider increasing X in"
                          " std::ratio<1, X> used for definition of model_step");
    throw std::invalid_argument(err);
  }

  if ((couplstep < condsubstep) |
      (couplstep < collsubstep))
  {
    const std::string err("sdm substepping warning: an sdm substep is"
                          " larger than the coupling step. "
                          " - it's possible but are you sure you want this?");
    throw std::invalid_argument(err);
  }

  if (couplstep < motionstep)
  {
    const std::string err("coupling step is smaller than the sdm motion step"
                          " - it's possible but are you sure you want this?");
    throw std::invalid_argument(err);
  }
}