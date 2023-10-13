/*
 * ----- CLEO -----
 * File: observers.hpp
 * Project: observers
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
 */


#ifndef OBSERVERS_HPP
#define OBSERVERS_HPP

#include <iostream>

#include "sdmdomain/gridbox.hpp"

struct Observer
{
private:
  unsigned int interval;

public:
  Observer(){}
  
  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void observe(const unsigned int t_mdl,
               Gridboxes &gbxs) const
  {
    std::cout << "obs gbxs @ t = " << t_mdl << "\n";
  }

  void observe_startstep(const unsigned int t_mdl,
                         Gridboxes &gbxs) const
  {
    if (on_step(t_mdl))
    {
      observe(t_mdl, gbxs);
    }
  }

  unsigned int get_obsstep() const
  {
    return interval;
  }
};

#endif // OBSERVERS_HPP