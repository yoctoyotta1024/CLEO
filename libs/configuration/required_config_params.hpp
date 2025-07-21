/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: required_config_params.hpp
 * Project: configuration
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for members of Config struct which determine CLEO's required configuration
 * parameters read from a config file.
 */

#ifndef LIBS_CONFIGURATION_REQUIRED_CONFIG_PARAMS_HPP_
#define LIBS_CONFIGURATION_REQUIRED_CONFIG_PARAMS_HPP_

#include <yaml-cpp/yaml.h>

#include <filesystem>
#include <iostream>
#include <limits>
#include <string>

/**
 * @brief Struct storing required configuration parameters for CLEO
 *
 * Required means parameters have no default values and therefore must be set upon
 * construction.
 *
 */
struct RequiredConfigParams {
  /* read configuration file given by config_filename to set members of required configuration */
  explicit RequiredConfigParams(const std::filesystem::path config_filename);

  void print_params() const;

  struct InputFilesParams {
    std::filesystem::path constants_filename; /**< filename for values of physical constants */
    std::filesystem::path grid_filename;      /**< filename for initialisation of GbxMaps */
  } inputfiles;

  struct OutputDataParams {
    std::filesystem::path setup_filename; /**< filename to copy model setup to */
    std::filesystem::path zarrbasedir;    /**< name of base directory of zarr output */
    size_t maxchunk;                      /**< maximum number of elements in zarr array chunks */
  } outputdata;

  struct DomainParams {
    unsigned int nspacedims; /**< no. of spatial dimensions to model */
    size_t ngbxs;            /**< total number of Gbxs */
    size_t maxnsupers;       /**< maximum number of SDs */
  } domain;

  struct TimestepsParams {
    double CONDTSTEP;   /**< time between SD condensation [s] */
    double COLLTSTEP;   /**< time between SD collision [s] */
    double MOTIONTSTEP; /**< time between SDM motion [s] */
    double COUPLTSTEP;  /**< time between thermodynamic couplings [s] */
    double OBSTSTEP;    /**< time between SDM observations [s] */
    double T_END;       /**< time span of integration from 0s to T_END [s] */
  } timesteps;
};

#endif  // LIBS_CONFIGURATION_REQUIRED_CONFIG_PARAMS_HPP_
