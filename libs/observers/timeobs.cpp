/*
 * ----- CLEO -----
 * File: timeobs.cpp
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
 * Functionality for bserver to output time
 * to array in zarr storage
 */


#include "./timeobs.hpp"

void TimeObs::at_start_step(const unsigned int t_mdl,
                    const viewh_constgbx h_gbxs) const
{
  std::cout << "obs time\n";
}