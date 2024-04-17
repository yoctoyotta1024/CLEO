/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: configparams.hpp
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

#ifndef LIBS_INITIALISE_CONFIGPARAMS_HPP_
#define LIBS_INITIALISE_CONFIGPARAMS_HPP_

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

  /*** SDM Runtime parameters ***/
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

/**
 * @brief Struct storing optional configuration parameters for CLEO
 *
 * Optional means parameters have default values and therefore need not be set upon
 * construction. Default values are not intended to be used and may caused model errors at runtime.
 *
 */
struct OptionalConfigParams {
  struct DoCondensationParams {
    using dblNaN = std::numeric_limits<double>::signaling_NaN;
    using uintNaN = std::numeric_limits<unsigned int>::signaling_NaN;
    bool do_alter_thermo = false; /**< enable condensation to alter the thermodynamic state */
    unsigned int iters = uintNaN; /**< suggested no. iterations of Newton Raphson Method */
    double SUBTSTEP = dblNaN;     /**< smallest timestep in cases where substepping occurs [s] */
    double rtol = dblNaN;         /**< relative tolerance for implicit Euler integration */
    double atol = dblNaN;         /**< abolute tolerance for implicit Euler integration */
  } condensation;

  struct FromFileDynamicsParams {
    std::string press_filename = ""; /**< name of file for pressure data */
    std::string temp_filename = "";  /**< name of file for temperature data */
    std::string qvap_filename = "";  /**< name of file for vapour mixing ratio data */
    std::string qcond_filename = ""; /**< name of file for liquid mixing ratio data */
    std::string wvel_filename = "";  /**< name of file for vertical (z) velocity data */
    std::string uvel_filename = "";  /**< name of file for horizontal x velocity data */
    std::string vvel_filename = "";  /**< name of file for horizontal y velocity data */
  } fromfiledynamics;

  struct CvodeDynamicsParams {
    using dblNaN = std::numeric_limits<double>::signaling_NaN;
    double P_INIT = dblNaN;    /**< initial pressure [Pa] */
    double TEMP_INIT = dblNaN; /**< initial temperature [T] */
    double relh_init = dblNaN; /**< initial relative humidity (%) */

    double W_AVG = dblNaN;  /**< average amplitude of w velocity sinusoid [m/s] (dP/dt ~ w*dP/dz) */
    double T_HALF = dblNaN; /**< timescale for w sinusoid, tau_half = T_HALF/pi [s] */
    double rtol = dblNaN;   /**< relative tolerance for integration of [P, T, qv, qc] ODEs */
    double atol = dblNaN;   /**< absolute tolerances for integration of [P, T, qv, qc] ODEs */
  } cvodedynamics;
};

#endif  // LIBS_INITIALISE_CONFIGPARAMS_HPP_
