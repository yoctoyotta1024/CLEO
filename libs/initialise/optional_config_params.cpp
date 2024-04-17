/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: optional_config_params.cpp
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
 * Functions involved in reading optional configuration parameters from a config file.
 */

#include "initialise/optional_config_params.hpp"

/* read configuration file given by config_filename to set members of required configuration */
OptionalConfigParams::OptionalConfigParams(const std::filesystem::path config_filename) {
  YAML::Node config = YAML::LoadFile(std::string{config_filename});

  // TODO(CB): assign optional values
}
