/*
 * ----- CLEO -----
 * File: timeobs.hpp
 * Project: observers
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 22nd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Observer to output gbxindx to array in a
 * zarr file system storage
 */

#ifndef TIMEOBS_HPP
#define TIMEOBS_HPP

#include <concepts>
#include <memory>
#include <functional>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/coordstorage.hpp"

namespace dlc = dimless_constants;

inline Observer auto
TimeObserver(const unsigned int interval,
             FSStore &store,
             const int maxchunk,
             const std::function<double(unsigned int)> step2dimlesstime);
/* constructs observer of time with a
constant timestep 'interval' using an
instance of the DoTimeObs class */

class DoTimeObs
/* observe time of 0th gridbox and write it
to an array 'zarr' store as determined by
the CoordStorage instance */
{
private:
  using store_type = CoordStorage<double>;
  std::shared_ptr<store_type> zarr;
  std::function<double(unsigned int)> step2dimlesstime; // function to convert timesteps to real time

public:
  DoTimeObs(FSStore &store,
          const int maxchunk,
          const std::function<double(unsigned int)> step2dimlesstime)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "time", "<f8",
                                          "s", dlc::TIME0)),
        step2dimlesstime(step2dimlesstime)
  {
    zarr->is_name("time");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes TimeObserver\n";
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  /* converts integer model timestep to dimensionless time,
  then writes to zarr coordinate storage */
  {
    const double time(step2dimlesstime(t_mdl));
    zarr->value_to_storage(time);
  }
};

inline Observer auto
TimeObserver(const unsigned int interval,
             FSStore &store,
             const int maxchunk,
             const std::function<double(unsigned int)> step2dimlesstime)
/* constructs observer of time with a
constant timestep 'interval' using an
instance of the DoTimeObs class */
{
  const auto obs = DoTimeObs(store, maxchunk, step2dimlesstime);
  return ConstTstepObserver(interval, obs);
}

#endif // TIMEOBS_HPP