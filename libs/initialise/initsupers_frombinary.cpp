/*
 * ----- CLEO -----
 * File: initsupers_frombinary.cpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 30th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for superdroplets' initial conditions
 * for CLEO SDM (e.g. superdroplet attributes)
 * by reading binary file. InitSupersFromBinary 
 * instance can be used by InitConds
 * struct as SuperdropInitConds type
 */


#include "./initsupers_frombinary.hpp"

InitSupersFromBinary::DataFromBinary
InitSupersFromBinary::fetch_data() const
{
  DataFromBinary initdata;

  initdata.solutes.at(0) = SoluteProperties{};


  return initdata;
}
