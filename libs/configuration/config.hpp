/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: config.hpp
 * Project: configuration
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for configuration Class including functions involved
 * in reading values from config files
 */

#ifndef LIBS_CONFIGURATION_CONFIG_HPP_
#define LIBS_CONFIGURATION_CONFIG_HPP_

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "configuration/copyfiles2txt.hpp"
#include "configuration/optional_config_params.hpp"
#include "configuration/required_config_params.hpp"

/**
 * @brief Struct storing configuration parameters read from a YAML file.
 *
 * This struct represents configuration settings for CLEO read in from
 * a configuration YAML file.
 */
struct Config {
 private:
  RequiredConfigParams required; /**< required configuration parameters of CLEO */
  OptionalConfigParams optional; /**< optional configuration parameters of CLEO */

 public:
  /**
   * @brief Constructor for Config.
   *
   * Initializes a Config instance by loading the configuration from the specified YAML
   * configuration file "config_filename". Then copy the setup to an output file "setup_filename".
   *
   * @param config_filename The name of the YAML configuration file.
   */
  explicit Config(const std::filesystem::path config_filename)
      : required(config_filename), optional(config_filename) {
    std::cout << "\n--- configuration ---\n";

    /* copy setup (config and constants files) to a txt file */
    const auto files2copy =
        std::vector<std::filesystem::path>{config_filename, required.inputfiles.constants_filename};
    copyfiles2txt(required.outputdata.setup_filename, files2copy);

    std::cout << "--- configuration: success ---\n";
  }

  std::filesystem::path get_grid_filename() const { return required.inputfiles.grid_filename; }

  std::filesystem::path get_zarrbasedir() const { return required.outputdata.zarrbasedir; }

  size_t get_maxchunk() const { return required.outputdata.maxchunk; }

  size_t get_maxnsupers() const { return required.domain.maxnsupers; }

  unsigned int get_nspacedims() const { return required.domain.nspacedims; }

  size_t get_ngbxs() const { return required.domain.ngbxs; }

  RequiredConfigParams::TimestepsParams get_timesteps() const { return required.timesteps; }

  Kokkos::InitializationSettings get_kokkos_initialization_settings() const {
    return optional.kokkos_settings.kokkos_initialization_settings;
  };

  OptionalConfigParams::CondensationParams get_condensation() const {
    return optional.condensation;
  }

  OptionalConfigParams::BreakupParams get_breakup() const { return optional.breakup; }

  OptionalConfigParams::InitSupersFromBinaryParams get_initsupersfrombinary() const {
    return optional.initsupersfrombinary;
  }

  OptionalConfigParams::CvodeDynamicsParams get_cvodedynamics() const {
    return optional.cvodedynamics;
  }

  OptionalConfigParams::FromFileDynamicsParams get_fromfiledynamics() const {
    return optional.fromfiledynamics;
  }

  OptionalConfigParams::YacDynamicsParams get_yac_dynamics() const { return optional.yac_dynamics; }

  OptionalConfigParams::AddSupersAtDomainTopParams get_addsupersatdomaintop() const {
    return optional.addsupersatdomaintop;
  }

  OptionalConfigParams::PythonBindingsParams get_python_bindings() const {
    return optional.python_bindings;
  }
};

#endif  // LIBS_CONFIGURATION_CONFIG_HPP_
