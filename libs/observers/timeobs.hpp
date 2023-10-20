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

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

class TimeObs
/* observe time of 0th gridbox and write it
to an array 'zarr' store as determined by
the CoordinateStorage instance */
{
private:
public:
  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const;
};

#endif // TIMEOBS_HPP