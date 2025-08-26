/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: optional_config_params.hpp
 * Project: configuration
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for members of Config struct which determine CLEO's optional configuration
 * parameters read from a config file.
 */

#ifndef LIBS_CONFIGURATION_OPTIONAL_CONFIG_PARAMS_HPP_
#define LIBS_CONFIGURATION_OPTIONAL_CONFIG_PARAMS_HPP_

#include <yaml-cpp/yaml.h>

#include <Kokkos_Core.hpp>
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

  void set_kokkos_settings(const YAML::Node& config);

  void print_kokkos_settings() const;

  void set_initsupers(const YAML::Node& config);

  void set_microphysics(const YAML::Node& config);

  void set_coupled_dynamics(const YAML::Node& config);

  void set_boundary_conditions(const YAML::Node& config);

  void set_python_bindings(const YAML::Node& config);

  /*** Kokkos Initialization Parameters ***/
  struct KokkosSettings {
    bool is_default = true; /**< true = default kokkos initialization */
    Kokkos::InitializationSettings kokkos_initialization_settings; /**< is default unless config */
  } kokkos_settings;

  /*** Super-Droplet Microphysics Parameters ***/
  struct CondensationParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    bool do_alter_thermo = false;        /**< true = cond/evap alters the thermodynamic state */
    size_t maxniters = NaNVals::sizet(); /**< maximum no. iterations of Newton Raphson Method */
    double MINSUBTSTEP = NaNVals::dbl(); /**< minimum subtimestep in cases of substepping [s] */
    double rtol = NaNVals::dbl();        /**< relative tolerance for implicit Euler integration */
    double atol = NaNVals::dbl();        /**< absolute tolerance for implicit Euler integration */
  } condensation;

  struct BreakupParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    struct ConstNFragsParams {
      double nfrags = NaNVals::dbl(); /**< average no. of fragments per droplet breakup */
    } constnfrags;
  } breakup;

  /*** Super-Droplet Initialisation Parameters ***/
  struct InitSupersFromBinaryParams {
    using fspath = std::filesystem::path;
    void set_params(const YAML::Node& config);
    void print_params() const;
    size_t maxnsupers = NaNVals::sizet();      /**< maximum number of SDs */
    fspath initsupers_filename = fspath();     /**< filename for initialisation of super-droplets */
    unsigned int nspacedims = NaNVals::uint(); /**< no. of spatial dimensions to model */
    size_t initnsupers = NaNVals::sizet();     /**< initial no. of super-droplets to initialise */
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

  struct YacDynamicsParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    double lower_longitude = NaNVals::dbl();
    double upper_longitude = NaNVals::dbl();
    double lower_latitude = NaNVals::dbl();
    double upper_latitude = NaNVals::dbl();
  } yac_dynamics;

  /*** Bounday Conditions Parameters ***/
  struct AddSupersAtDomainTopParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    size_t initnsupers = NaNVals::sizet(); /**< initial no. of super-droplets in domain */
    size_t newnsupers = NaNVals::sizet();  /**< number SDs to add to each gridbox above COORD3LIM */
    double COORD3LIM = NaNVals::dbl();     /**< SDs added to domain with coord3 >= COORD3LIM [m] */
    double DRYRADIUS = NaNVals::dbl();     /**< dry radius of new super-droplets (for msol) [m] */
    double MINRADIUS = NaNVals::dbl();     /**< minimum radius of new super-droplets [m] */
    double MAXRADIUS = NaNVals::dbl();     /**< maximum radius of new super-droplets [m] */
    double NUMCONC_a = NaNVals::dbl();     /**< number conc. of 1st droplet lognormal dist [m^-3] */
    double GEOMEAN_a = NaNVals::dbl();     /**< geometric mean radius of 1st lognormal dist [m] */
    double geosigma_a = NaNVals::dbl(); /**< geometric standard deviation of 1st lognormal dist */
    double NUMCONC_b = NaNVals::dbl();  /**< number conc. of 2nd droplet lognormal dist [m^-3] */
    double GEOMEAN_b = NaNVals::dbl();  /**< geometric mean radius of 2nd lognormal dist [m] */
    double geosigma_b = NaNVals::dbl(); /**< geometric standard deviation of 2nd lognormal dist */
  } addsupersatdomaintop;

  /** CLEO Python Bindings Parameters */
  struct PythonBindingsParams {
    void set_params(const YAML::Node& config);
    void print_params() const;
    bool enable_terminal_velocity =
        false;                        /**< true enables terminal velocity in superdroplet motion */
    bool enable_condensation = false; /**< true enables condensation in microphysics */
    bool enable_collisions = false;   /**< true enables collisions in microphysics */
    struct Observers {
      bool time = false;
      bool gbxindex = true;
      bool totnsupers = false;
      bool massmoms = false;
      bool rainmassmoms = false;
      bool gridboxes = false;
      bool superdrops = false;
      bool precip = false;
    } enable_observers; /**< true for set of booleans in struct enables various observers */
  } python_bindings;
};

#endif  // LIBS_CONFIGURATION_OPTIONAL_CONFIG_PARAMS_HPP_
