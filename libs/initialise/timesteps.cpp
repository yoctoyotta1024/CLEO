/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: timesteps.cpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * File Description:
 * functionality for handling model timesteps and
 * their conversions to/from real times
 */

#include "initialise/timesteps.hpp"

/* (dimensionless) double's that are timesteps in config struct
are converted into integer values of model timesteps using
model_step and secd template functions created using std::chrono library.
Throw error if after convertion into model timestep, any
timestep = 0. Substeps for sdmprocess must be larger than steps! */
Timesteps::Timesteps(const RequiredConfigParams::TimestepsParams &config)
    : condstep(realtime2step(config.CONDTSTEP)),
      collstep(realtime2step(config.COLLTSTEP)),
      motionstep(realtime2step(config.MOTIONTSTEP)),
      couplstep(realtime2step(config.COUPLTSTEP)),
      obsstep(realtime2step(config.OBSTSTEP)),
      t_end(realtime2step(config.T_END)) {
  if ((condstep == 0) | (collstep == 0) | (motionstep == 0) | (couplstep == 0) | (obsstep == 0) |
      (t_end == 0)) {
    const std::string err(
        "A model step = 0, possibly due to bad"
        "conversion of a real timestep [s]. Consider"
        " increasing X in std::ratio<1, X> used for"
        " definition of model_step");
    throw std::invalid_argument(err);
  }

  const auto maxsubstep = std::max(condstep, collstep);
  const auto minstep = std::min(std::min(couplstep, obsstep), motionstep);
  if (minstep < maxsubstep) {
    const std::string err(
        "invalid microphysics sub-stepping: "
        "a microphysics substep is greater "
        "motion / coupling / observation step");
    throw std::invalid_argument(err);
  }

  if (std::min(couplstep, obsstep) < motionstep) {
    const std::string err(
        "invalid CLEO SDM sub-stepping: "
        "motion / microphysics step is greater "
        "than coupling / observation step");
    throw std::invalid_argument(err);
  }
}
