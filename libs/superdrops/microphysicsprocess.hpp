/*
 * ----- CLEO -----
 * File: microphysicsprocess.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 14th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Microphysical Process Concept as well as structures for
 * various processes that may occur in CLEO SDM,
 * eg. condensation or collision-coalescence
 * (see ConstTstepProcess struct)
 */

#ifndef MICROPHYSICSPROCESS_HPP
#define MICROPHYSICSPROCESS_HPP

#include <iostream>

struct MicrophysicsProcess
{
  MicrophysicsProcess(){}

  unsigned int next_step(const unsigned int t_sdm)const
  {
    return t_sdm + 100;
  }
  
  void run_step(const unsigned int subt) const
  {
    std::cout << "microphys @ t = " << subt << "\n";
  }
};

#endif // MICROPHYSICSPROCESS_HPP