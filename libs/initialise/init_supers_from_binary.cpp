/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: init_supers_from_binary.cpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 19th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct for reading in some super-droplets' initial conditions for CLEO SDM
 * (e.g. superdroplet attributes) from a binary file. InitAllSupersFromBinary instance
 * can be used by InitConds struct as SuperdropInitConds type.
 */

#include "initialise/init_supers_from_binary.hpp"

InitSupersData InitSupersFromBinary::fetch_invalid_superdrops_data(InitSupersData &initdata) const {
  std::cout << "size: " << maxnsupers << ", " << initdata.radii.size() << "\n";
  return initdata;
}
