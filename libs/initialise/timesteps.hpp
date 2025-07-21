/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: timesteps.hpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Structs for handling model timesteps and
 * their conversions to/from real times
 */

#ifndef LIBS_INITIALISE_TIMESTEPS_HPP_
#define LIBS_INITIALISE_TIMESTEPS_HPP_

#include <algorithm>
#include <chrono>
#include <ratio>
#include <stdexcept>
#include <string>

#include "../cleoconstants.hpp"
#include "configuration/required_config_params.hpp"

namespace dlc = dimless_constants;

using model_step = std::chrono::duration<int, std::ratio<1, 100000>>;
using secd = std::chrono::duration<double>;

/* convert TSTEP [seconds] (e.g. a double from Config struct)
into a dimensionless time and then into an integer no. of
model steps using model_step chrono function */
inline unsigned int realtime2step(const double TSTEP) {
  auto step(std::chrono::round<model_step>(secd{TSTEP / dlc::TIME0}));
  return step.count();
}

/* convert TSTEP [seconds] (e.g. a double from Config
struct) into a dimensionless time  */
inline double realtime2dimless(const double TSTEP) { return TSTEP / dlc::TIME0; }

/* convert model step (integer) into a time [seconds]
given secd and model_step chrono functions */
inline double step2realtime(const unsigned int mdlstep) {
  return secd{model_step{mdlstep}}.count() * dlc::TIME0;
}

/* convert model timestep (integer) into a dimensionless
time given secd and model_step chrono functions */
inline double step2dimlesstime(const unsigned int mdlstep) {
  return secd{model_step{mdlstep}}.count();
}

/* integer intervals (timesteps) involved in running CLEO model */
class Timesteps {
 private:
  const unsigned int condstep;
  const unsigned int collstep;
  const unsigned int motionstep;
  const unsigned int couplstep;
  const unsigned int obsstep;
  const unsigned int t_end;

 public:
  /* (dimensionless) double's that are timesteps in config struct
  are converted into integer values of model timesteps using
  model_step and secd template functions created using
  std::chrono library. Throw error if after convertion into
  model timestep, any timestep = 0 or if a sub-timestep is
  longer than a timestep */
  explicit Timesteps(const RequiredConfigParams::TimestepsParams &config_tsteps);

  auto get_condstep() const { return condstep; }
  auto get_collstep() const { return collstep; }
  auto get_motionstep() const { return motionstep; }
  auto get_couplstep() const { return couplstep; }
  auto get_obsstep() const { return obsstep; }
  auto get_t_end() const { return t_end; }
};

#endif  // LIBS_INITIALISE_TIMESTEPS_HPP_
