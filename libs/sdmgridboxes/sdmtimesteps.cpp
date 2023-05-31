// Author: Clara Bayley
// File: sdmtimesteps.cpp
/* Structs for handling values of
timstep variables for SDM */

#include "sdmtimesteps.hpp"

SDMTimesteps::SDMTimesteps(const double CONDTSTEP, const double COLLTSTEP,
                           const double MOTIONTSTEP, const double COUPLTSTEP,
                           const double OBSTSTEP, const double T_END)
    /* (dimensionless) double's that are timesteps in config struct
    are converted into integer values of model timesteps using
    model_step and secd template functions created using std::chrono library.
    Throw error if after convertion into model timestep, any
    timestep = 0. Substeps for sdmprocess must be larger than steps! */
    : condsubstep(realtime2step(CONDTSTEP)),
      collsubstep(realtime2step(COLLTSTEP)),
      motionstep(realtime2step(MOTIONTSTEP)),
      couplstep(realtime2step(COUPLTSTEP)),
      obsstep(realtime2step(OBSTSTEP)),
      t_end(realtime2step(T_END))
{
  if ((condsubstep == 0) | (collsubstep == 0) | (motionstep == 0) |
      (couplstep == 0) | (obsstep == 0) | (t_end == 0))
  {
    const std::string err("A model step = 0, possibly due to bad"
                          "conversion of a real timestep [s]. Consider"
                          " increasing X in std::ratio<1, X> used for"
                          " definition of model_step");
    throw std::invalid_argument(err);
  }

  const int minstep = std::min(std::min(couplstep, obsstep), motionstep);
  if ((minstep < condsubstep) |
      (minstep < collsubstep))
  {
    const std::string err("invalid sdm substepping: an sdm substep"
                          " is larger than the smallest step"
                          " (coupling, observation or motion step)");
    throw std::invalid_argument(err);
  }

  if (std::min(couplstep, obsstep) < motionstep)
  {
    const std::string err("Warning: coupling / observation step is"
                          "smaller than the sdmmotion step"
                          " - are you really sure you want this?");
    throw std::invalid_argument(err);
  }
}