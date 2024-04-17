/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: config.hpp
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
 * Header file for configuration Class including functions involved
 * in reading values from config files
 */

#ifndef LIBS_INITIALISE_CONFIG_HPP_
#define LIBS_INITIALISE_CONFIG_HPP_

#include <yaml-cpp/yaml.h>

#include <string_view>

#include "./copyfiles2txt.hpp"

/**
 * @brief Struct storing configuration parameters read from a YAML file.
 *
 * This struct represents configuration settings for CLEO read in from
 * a configuration YAML file.
 */
struct Config {
 private:
 public:
  /**
   * @brief Constructor for Config.
   *
   * Initializes a Config instance by loading the configuration
   * from the specified YAML configuration file and copies the setup
   * to an output file "setup_filename".
   *
   * @param config_filename The name of the YAML configuration file.
   */
  explicit Config(const std::string_view config_filename) {
    std::cout << "\n--- configuration ---\n";

    loadconfiguration(config_filename);

    /* copy setup (config and constants files) to a txt file */
    const std::string filestr(config_filename);
    copyfiles2txt(setup_filename, {filestr, constants_filename});

    std::cout << "--- configuration: success ---\n";
  }
};

#endif  // LIBS_INITIALISE_CONFIG_HPP_
