/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: required_config_params.hpp
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
 * Header file for members of Config struct which determine CLEO's configuration
 */

#ifndef LIBS_INITIALISE_REQUIRED_CONFIG_PARAMS_HPP_
#define LIBS_INITIALISE_REQUIRED_CONFIG_PARAMS_HPP_

#include <filesystem>
#include <limits>
#include <string>
#include <string_view>

// TODO(CB): check types of config params e.g. int maxchunk -> size_t maxchunk
// TODO(CB): use std::filesystem::path not string

/**
 * @brief Struct storing required configuration parameters for CLEO
 *
 * Required means parameters have no default values and therefore must be set upon
 * construction.
 *
 */
struct RequiredConfigParams {
  /*** Input and Output Parameters ***/
  struct InputFilesParams {
    std::string constants_filename;  /**< name of input file for values of physical constants */
    std::string initsupers_filename; /**< name of input file for initialisation of super-droplets */
    std::string grid_filename;       /**< name of input file for initialisation of GbxMaps */
  } inputfiles;

  struct OutputDataParams {
    std::string setup_filename;        /**< name of output file to copy model setup to */
    std::string stats_filename;        /**< name of output file to output runtime statistics to */
    std::filesystem::path zarrbasedir; /**< name of base directory of zarr output */
    size_t maxchunk;                   /**< maximum number of elements in zarr array chunks */
  } outputdata;

  /*** SDM Runtime Parameters ***/
  struct DomainParams {
    unsigned int nspacedims;           /**< no. of spatial dimensions to model */
    size_t ngbxs;                      /**< total number of Gbxs */
    size_t totnsupers;                 /**< initial total no. of SDs */
    std::string_view coupled_dynamics; /**< type of coupled dynamics to configure */
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

#endif  // LIBS_INITIALISE_REQUIRED_CONFIG_PARAMS_HPP_
