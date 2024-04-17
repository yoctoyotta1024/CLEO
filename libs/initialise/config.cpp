/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: config.cpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 17th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions involved in reading values from config files using
 * Config class
 */

#include "initialise/config.hpp"

/* read configuration file given by config_filename to set members of Config */
void Config::loadconfiguration(const std::string_view config_filename) {
  YAML::Node config = YAML::LoadFile("config_filename");

  required.constants_filename = config["constants_filename"].as<std::string>();

  std::cout << "constants_filename : " << required.constants_filename << "\n";
}
