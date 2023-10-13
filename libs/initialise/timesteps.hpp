/*
 * ----- CLEO -----
 * File: timesteps.hpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Structs for handling model timesteps and
 * their conversions to/from real times
 */

#ifndef TIMESTEPS_HPP
#define TIMESTEPS_HPP

#include <chrono>
#include <ratio>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "./config.hpp"
#include "../cleoconstants.hpp"

namespace dlc = dimless_constants;

using model_step = std::chrono::duration<int, std::ratio<1, 100000>>;
using secd = std::chrono::duration<double>;

inline unsigned int realtime2step(const double TSTEP)
/* convert TSTEP [seconds] (e.g. a double from Config struct)
into a dimensionless time and then into an integer no. of
model steps using model_step chrono function */
{
  auto step(std::chrono::round<model_step>(secd{TSTEP / dlc::TIME0}));
  return step.count();
}

inline double realtime2dimless(const double TSTEP)
/* convert TSTEP [seconds] (e.g. a double from Config
struct) into a dimensionless time  */
{
  return TSTEP / dlc::TIME0;
}

inline double step2realtime(const unsigned int mdlstep)
/* convert model step (integer) into a time [seconds]
given secd and model_step chrono functions */
{
  return secd{model_step{mdlstep}}.count() * dlc::TIME0;
}

inline double step2dimlesstime(const unsigned int mdlstep)
/* convert model timestep (integer) into a dimensionless
time given secd and model_step chrono functions */
{
  return secd{model_step{mdlstep}}.count();
}

class Timesteps
/* integer intervals (timesteps) involved in running CLEO model */
{
private:
  const unsigned int condstep;
  const unsigned int collstep;
  const unsigned int motionstep;
  const unsigned int couplstep;
  const unsigned int obsstep;
  const unsigned int t_end;

public:
  Timesteps(const Config &config);
  /* (dimensionless) double's that are timesteps in config struct
  are converted into integer values of model timesteps using
  model_step and secd template functions created using
  std::chrono library. Throw error if after convertion into
  model timestep, any timestep = 0 or if a sub-timestep is
  longer than a timestep */

  unsigned int get_condstep() const { return condstep; }
  unsigned int get_collstep() const { return collstep; }
  unsigned int get_motionstep() const { return motionstep; }
  unsigned int get_couplstep() const { return couplstep; }
  unsigned int get_obsstep() const { return obsstep; }
  unsigned int get_t_end() const { return t_end; }
};

#endif // TIMESTEPS_HPP