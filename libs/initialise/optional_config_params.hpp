/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: optional_config_params.hpp
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

#ifndef LIBS_INITIALISE_OPTIONAL_CONFIG_PARAMS_HPP_
#define LIBS_INITIALISE_OPTIONAL_CONFIG_PARAMS_HPP_

#include <filesystem>
#include <limits>
#include <string>
#include <string_view>

// TODO(CB): check types of config params e.g. int maxchunk -> size_t maxchunk
// TODO(CB): use std::filesystem::path not string

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

#endif  // LIBS_INITIALISE_OPTIONAL_CONFIG_PARAMS_HPP_
