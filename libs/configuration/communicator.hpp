/*
 * Copyright (c) 2025 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: communicator.hpp
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
 * Header file for members of Config struct which determine CLEO's required configuration
 * parameters read from a config file.
 */

#ifndef LIBS_CONFIGURATION_COMMUNICATOR_HPP_
#define LIBS_CONFIGURATION_COMMUNICATOR_HPP_

#include <mpi.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "configuration/config.hpp"

class init_communicator {
 private:
  bool yac_present;
  // const Config &config;
 public:
  int rank;
  int size;
  MPI_Comm comm;
  // Constructor
  // init_communicator(const std::filesystem::path config_filename)
  explicit init_communicator(const Config &config);
};

#endif  // LIBS_CONFIGURATION_COMMUNICATOR_HPP_
