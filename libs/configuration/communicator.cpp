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
int init_communicator::comm_size = -1;
int init_communicator::my_rank = -1;

init_communicator::init_communicator(int argc, char* argv[], const Config& config) {
  is_using_yac = config.get_yac_settings().is_using_yac;
  if (is_using_yac) {
    std::cout << "yac is present\n";
    // -- YAC initialization and calendar definitions ---
    yac_cinit();
    yac_cread_config_yaml(config.get_yac_settings().yac_config_file.string().c_str());
    // --- Component definition ---
    const auto cleo_component_name = config.get_yac_settings().cleo_component_name;
    yac_cdef_comp(cleo_component_name.c_str(), &init_communicator::yac_comp_id);
    yac_cget_comp_comm(init_communicator::yac_comp_id, &comm);
    MPI_Comm_size(comm, &init_communicator::comm_size);
    MPI_Comm_rank(comm, &my_rank);
  } else {
    std::cout << "yac is not present " << is_using_yac << "\n";

    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
      MPI_Init(&argc, &argv);
      MPI_Initialized(&mpi_initialized);
    }
    std::cout << "MPI initialized " << mpi_initialized << "\n";

    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &my_rank);
  }
};

init_communicator::~init_communicator() {
  if (is_using_yac) {
    std::cout << "yac_finalized elsewhere" << std::endl;
  } else {
    std::cout << "mpi finalizing" << std::endl;
    MPI_Finalize();
  }
};

MPI_Comm init_communicator::get_communicator() {
  if (init_communicator::comm == MPI_COMM_NULL) {
    std::cout << "Communicator not initialized, calling MPI Abort!" << std::endl;
    MPI_Abort(comm, 1);
  }
  return comm;
};

int init_communicator::get_yac_comp_id() {
  if (!(yac_comp_id > 0)) {
    std::cout << "Invalid yac_comp_id, calling MPI Abort!" << std::endl;
    MPI_Abort(comm, 1);
  }
  return init_communicator::yac_comp_id;
};

int init_communicator::get_comm_size() { return init_communicator::comm_size; };

int init_communicator::get_comm_rank() { return my_rank; };
