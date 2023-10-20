/*
 * ----- CLEO -----
 * File: timeobs.hpp
 * Project: observers
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Observer to output time to array in zarr storage
 */

#ifndef TIMEOBS_HPP
#define TIMEOBS_HPP

#include <iostream>
#include <concepts>
#include <memory>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/coordstorage.hpp"

class TimeObs
/* observe time of 0th gridbox and write it
to an array 'zarr' store as determined by
the CoordStorage instance */
{
private:
  using store_type = CoordStorage<double>;
  std::shared_ptr<store_type> zarr;
  std::function<double(int)> step2dimlesstime; // function to convert timesteps to real time

public:
  TimeObs(FSStore &store,
          const int maxchunk,
          const std::function<double(int)> step2dimlesstime)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "time", "<f8",
                                          "s", dlc::TIME0)),
        step2dimlesstime(step2dimlesstime)
  {
    zarr->is_name("time");
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
             const std::function<double(int)> step2dimlesstime)
/* constructs Microphysical Process for
condensation/evaporation of superdroplets with a
constant timestep 'interval' given the
"do_condensation" function-like type */
{
  const auto obs(TimeObs(store, maxchunk, step2dimlesstime));
  return ConstTstepObserver(interval, obs);
}

#endif // TIMEOBS_HPP