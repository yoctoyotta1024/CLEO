// Author: Clara Bayley
// File: timesteps.cpp
/* Structs for handling values of
timstep variables for SDM */

#include "timesteps.hpp"


Timesteps::Timesteps(Config config)
/* (dimensionless) double's that are timesteps in config struct
are converted into integer values of model timesteps using
model_step and secd template functions created using std::chrono library.
Throw error if after convertion into model timestep, any
timestep = 0 */
    : condstep(realtime2timestep(config.COND_TSTEP)),
      collstep(realtime2timestep(config.COLL_TSTEP)),
      sedistep(realtime2timestep(config.SEDI_TSTEP)),
      xchangestep(realtime2timestep(config.XCHANGE_TSTEP)),
      outstep(realtime2timestep(config.OUT_TSTEP)),
      tend(realtime2timestep(config.TEND))
{
  if ((condstep == 0) | (collstep == 0) | (sedistep == 0) |
      (xchangestep == 0) | (outstep == 0) | (tend == 0))
  {
    throw std::invalid_argument("ERROR! A MODEL TIMESTEP = 0.\n"
                                "Consider increasing model_step's denom "
                                "in std::ratio<1, denom>>");
  }

  if ((xchangestep < condstep) |
      (xchangestep < collstep) |
      (xchangestep < sedistep))
  {
    throw std::invalid_argument("ERROR! XCHANGE MODEL TIMESTEP less than"
                                "cond, coll or sedi tstep. undefined sdm "
                                "timestepping (this makes no sense)");
  }

  if ((outstep < condstep) |
      (outstep < collstep) |
      (outstep < sedistep))
  {
    throw std::invalid_argument("ERROR! OUTSTEP MODEL TIMESTEP less than cond"
                                "coll or sedi tstep. undefined sdm timstepping");
  }

}