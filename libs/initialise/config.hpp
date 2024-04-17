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

#include "./configparams.hpp"
#include "./copyfiles2txt.hpp"

/**
 * @brief Struct storing configuration parameters read from a YAML file.
 *
 * This struct represents configuration settings for CLEO read in from
 * a configuration YAML file.
 */
struct Config {
 private:
  /* read configuration file given by config_filename to set members of Config */
  void loadconfiguration(const std::string_view config_filename);

  RequiredConfigParams required; /**< required configuration parameters of CLEO */

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
    const auto files2copy =
        std::vector<std::string>{std::string{config_filename}, required.constants_filename};
    copyfiles2txt(required.setup_filename, files2copy);

    std::cout << "--- configuration: success ---\n";
  }

  std::string get_initsupers_filename() const { return required.initsupers_filename; };
  std::string get_grid_filename() const { return required.grid_filename; }
  std::string get_stats_filename() const { return required.stats_filename; }
  std::filesystem::path get_zarrbasedir() const { return required.zarrbasedir; }
  size_t get_maxchunk() const { return required.maxchunk; }
  unsigned int get_nspacedims() const { return required.nspacedims; }
  size_t get_ngbxs() const { return required.ngbxs; }
  size_t get_totnsupers() const { return required.totnsupers; }
};

#endif  // LIBS_INITIALISE_CONFIG_HPP_
