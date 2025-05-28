/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: communicator.cpp
 * Project: configuration
 * Created Date: Tuesday 06 May 2025
 * Author: Lakshmi Aparna Devulapalli (LAD)
 * Additional Contributors: Clara Bayley (CB)
 * -----
 * Last Modified: Wednesday 28th May 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality for members of Config struct which determine CLEO's required configuration
 * parameters read from a config file.
 */

#include "configuration/communicator.hpp"

extern "C" {
#include "yac.h"
}

int init_communicator::yac_comp_id = -1;
MPI_Comm init_communicator::comm = NULL;

init_communicator::init_communicator(const Config &config) {
  if (!(std::isnan(config.get_yac_dynamics().lower_longitude))) {
    std::cout << "yac is present";
    // -- YAC initialization and calendar definitions ---
    yac_cinit();
    // --- Component definition ---
    std::string component_name = "cleo";
    yac_cdef_comp(component_name.c_str(), &init_communicator::yac_comp_id);
    yac_cget_comp_comm(init_communicator::yac_comp_id, &comm);
  } else {
    std::cout << "yac is not present" << yac_present;
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    int comm_size;
    MPI_Comm_size(comm, &comm_size);
    if (comm_size > 1) {
      std::cout << "ERROR: The current example is not prepared"
                << " to be run with more than one MPI process" << std::endl;
      MPI_Abort(comm, 1);
    }
  }
};

MPI_Comm init_communicator::get_communicator() {
  if (init_communicator::comm == MPI_COMM_NULL) {
    // call MPI_Abort and print msg
  }
  return comm;
};
