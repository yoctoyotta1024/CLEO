/*
 * ----- CLEO -----
 * File: condensation.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
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
 * struct for modelling condensation
 * microphysical process in SDM
 */

#ifndef CONDENSATION_HPP
#define CONDENSATION_HPP

#include <iostream>

struct Condensation
{
private:
  unsigned int interval;

  void do_condensation(const unsigned int subt) const
  {
    std::cout << "cond microphys @ t = " << subt << "\n";
  }

public:
  Condensation(const unsigned int condstep)
      : interval(condstep) {}

  unsigned int next_step(const unsigned int t_sdm) const
  {
    return ((t_sdm / interval) + 1) * interval;
  }
  
  bool on_step(const unsigned int t_sdm) const
  {
    return t_sdm % interval == 0;
  }

  void run_step(const unsigned int subt) const
  {
    do_condensation(subt);
  }
};

#endif // CONDENSATION_HPP