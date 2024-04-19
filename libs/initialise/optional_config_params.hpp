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
 * Last Modified: Friday 19th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for members of Config struct which determine CLEO's optional configuration
 * parameters read from a config file.
 */

#ifndef LIBS_INITIALISE_OPTIONAL_CONFIG_PARAMS_HPP_
#define LIBS_INITIALISE_OPTIONAL_CONFIG_PARAMS_HPP_

#include <yaml-cpp/yaml.h>

#include <cassert>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>

namespace NaNVals {
inline double dbl() { return std::numeric_limits<double>::signaling_NaN(); };
inline unsigned int uint() { return std::numeric_limits<unsigned int>::signaling_NaN(); };
inline size_t sizet() { return std::numeric_limits<size_t>::signaling_NaN(); };
}  // namespace NaNVals

/**
 * @brief Struct storing optional configuration parameters for CLEO
 *
 * Optional means parameters have default values and therefore need not be set upon
 * construction. Default values are not intended to be used and may caused model errors at runtime.
 *
 */
struct OptionalConfigParams {
  /* read configuration file given by config_filename to set members of required configuration */
  explicit OptionalConfigParams(const std::filesystem::path config_filename);

  void set_initsupers(const YAML::Node& config);

  void set_microphysics(const YAML::Node& config);

  void set_coupled_dynamics(const YAML::Node& config);

  void set_boundary_conditions(const YAML::Node& config);

  /*** Super-Droplet Microphysics Parameters ***/
  struct CondensationParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    bool do_alter_thermo = false;          /**< true = cond/evap alters the thermodynamic state */
    unsigned int niters = NaNVals::uint(); /**< suggested no. iterations of Newton Raphson Method */
    double SUBTSTEP = NaNVals::dbl();      /**< smallest subtimestep in cases of substepping [s] */
    double rtol = NaNVals::dbl();          /**< relative tolerance for implicit Euler integration */
    double atol = NaNVals::dbl();          /**< abolute tolerance for implicit Euler integration */
  } condensation;

  /*** Super-Droplet Initialisation Parameters ***/
  struct InitSupersFromBinaryParams {
    using fspath = std::filesystem::path;
    void set_params(const YAML::Node& config);
    void print_params() const;
    fspath initsupers_filename = fspath();     /**< filename for initialisation of super-droplets */
    unsigned int nspacedims = NaNVals::uint(); /**< no. of spatial dimensions to model */
    size_t initnsupers = NaNVals::sizet();     /**< no. of super-droplets to initialise */
  } initsupersfrombinary;

  /*** Coupled Dynamics Parameters ***/
  struct FromFileDynamicsParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    using fspath = std::filesystem::path;
    unsigned int nspacedims = NaNVals::uint(); /**< no. of spatial dimensions to model */
    fspath press = fspath();                   /**< name of file for pressure data */
    fspath temp = fspath();                    /**< name of file for temperature data */
    fspath qvap = fspath();                    /**< name of file for vapour mixing ratio data */
    fspath qcond = fspath();                   /**< name of file for liquid mixing ratio data */
    fspath wvel = fspath();                    /**< name of file for vertical (z) velocity data */
    fspath uvel = fspath();                    /**< name of file for horizontal x velocity data */
    fspath vvel = fspath();                    /**< name of file for horizontal y velocity data */
  } fromfiledynamics;

  struct CvodeDynamicsParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    size_t ngbxs = NaNVals::sizet(); /**< no. of spatial dimensions to model */
    /* initial (uniform) thermodynamic conditions */
    double P_init = NaNVals::dbl();    /**< initial pressure [Pa] */
    double TEMP_init = NaNVals::dbl(); /**< initial temperature [T] */
    double relh_init = NaNVals::dbl(); /**< initial relative humidity (%) */
    /* ODE solver parameters */
    double W_avg = NaNVals::dbl(); /**< average amplitude of sinusoidal w [m/s] (dP/dt ~ w*dP/dz) */
    double TAU_half = NaNVals::dbl(); /**< timescale for w sinusoid, tau_half = TAU_half/pi [s] */
    double rtol = NaNVals::dbl(); /**< relative tolerance for integration of [P, T, qv, qc] ODEs */
    double atol = NaNVals::dbl(); /**< absolute tolerances for integration of [P, T, qv, qc] ODEs */
  } cvodedynamics;

  /*** Bounday Conditions Parameters ***/
  struct AddSupersAtDomainTopParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    double COORD3LIM = NaNVals::dbl();    /**< SDs added to domain with coord3 >= COORD3LIM [m] */
    size_t newnsupers = NaNVals::sizet(); /**< number SDs to add to each gridbox above COORD3LIM */
  } addsupersatdomaintop;
};

#endif  // LIBS_INITIALISE_OPTIONAL_CONFIG_PARAMS_HPP_
