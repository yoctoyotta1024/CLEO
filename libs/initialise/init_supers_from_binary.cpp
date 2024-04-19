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

/* returns InitSupersData created by reading some data from a binary file and
filling the rest with invalid super-droplets */
InitSupersData InitSupersFromBinary::fetch_superdrops_from_file() const {
  auto initdata = InitSupersData{};  // TODO(CB) WIP

  return initdata
}

InitSupersData InitSupersFromBinary::fetch_invalid_superdrops_data(InitSupersData &initdata) const {
  // tTODO(CB): WIP
  return initdata;
}
