/*
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: superdrops
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 16th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates given
 * a formula for its terminal velocity and the wind
 * velocity obtained via a simple linear interpolation.
 * Methods follows equations in Grabowski et al. 2018
 */


#ifndef PREDCORRMOTION_HPP
#define PREDCORRMOTION_HPP

#include <iostream>

struct PredCorrMotion
{
private:
  unsigned int interval;

public:
  PredCorrMotion(const unsigned int motionstep)
      : interval(motionstep) {}

  unsigned int next_step(const unsigned int t_sdm) const
  {
    return ((t_sdm / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_sdm) const
  {
    return t_sdm % interval == 0;
  }

  void update_superdrop_coords(const unsigned int t_sdm) const;

};

#endif // PREDCORRMOTION_HPP