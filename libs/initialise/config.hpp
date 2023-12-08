/*
 * ----- CLEO -----
 * File: config.hpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 8th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Header file for configuration Class including functions involved
 * in reading values from config files
 */

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <fstream>
#include <string>
#include <string_view>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <limits>

#include "./copyfiles2txt.hpp"

struct Config
{
private:
  struct NameValue
  {
    std::string name;
    std::string value;
  };

  void loadconfiguration(const std::string_view config_filename);
  /* read configuration file given by config_filename
  then calls configvariable to assign value to corresponding
  member of Config structure */

  void open_file(const std::string_view config_filename,
                 std::ifstream &file) const;
  /* opens file or throws error if opening fails */

  NameValue getvariable_fromline(const std::string line) const;
  /* extracts name of variable and it's value from a string
  of the form 'variable_name = value    #comment' */

  bool string2bool(const std::string value) const;
  /* if string == 'true', 'True' or '1', returns true
  C++ boolean, else returns false */

  void configvariable(const std::string name, const std::string value);
  /* setter function. assigns value of member of Config struct
  called 'name' by coverting strings containing
  value into actual value for that  members's type */

  void configvariable_thermosolverfromfile(const std::string name,
                                           const std::string value);
  /* setter function for assigning 'value' to members of
  Config struct called 'name' specifically for members
  involved when thermosolver == 'fromfile' */

  void configvariable_thermosolvercvode(const std::string name,
                                        const std::string value);
  /* setter function for assigning 'value' to members of
  Config struct called 'name' specifically for members
  involved when thermosolver == 'cvode' */

public:
  /*** Initialisation and Output Data parameters ***/
  std::string constants_filename;    // filename containing values of physical constants
  std::string initsupers_filename;   // binary filename for initialisation of SDs
  std::string grid_filename;         // binary filename for GBx boundaries
  std::string setuptxt;              // name of .txt output file to copy setup to
  std::string stats_filename;        // name of .txt file to output runtime statistics to
  std::filesystem::path zarrbasedir; // zarr store base directory
  int maxchunk;                      // maximum no. of elements in chunks of zarr store array

  /*** SDM Runtime parameters ***/
  /* domain setup */
  int nspacedims; // no. of spatial dimensions to model
  int ngbxs;      // total number of Gbxs
  int totnsupers; //(initial) total no. of SDs

  /* timestepping */
  double CONDTSTEP;   // time between SD condensation events [s]
  double COLLTSTEP;   // time between SD collision events [s]
  double MOTIONTSTEP; // time between SDM motion [s]
  double COUPLTSTEP;  // time between thermodynamic couplings [s]
  double OBSTSTEP;    // time between SDM observations [s]
  double T_END;       // time span of integration from 0s to T_END [s]

  /* microphysics */
  unsigned int cond_iters; // suggested no. iterations of Newton Raphson Method
  double cond_SUBTSTEP;    // smallest timestep in cases where substepping occurs [s]
  double cond_rtol;        // relative tolerance for implicit euler integration
  double cond_atol;        // abolute tolerance for implicit euler integration

  /* superdroplets */
  bool doAlterThermo; // enable condensation to alter the thermodynamic state

  /*** Coupled Dynamics Solver Parameters ***/
  std::string thermosolver; // type of dynamics solver to configure

  /* read in dynamics from file (default to empty) */
  std::string press_filename = ""; // binary filename for pressure
  std::string temp_filename = "";  // binary filename for temperature
  std::string qvap_filename = "";  // binary filename for vapour mixing ratio
  std::string qcond_filename = ""; // binary filename for liquid mixing ratio
  std::string wvel_filename = "";  // binary filename for vertical (z) velocity
  std::string uvel_filename = "";  // binary filename for horizontal x velocity
  std::string vvel_filename = "";  // binary filename for horizontal y velocity

  /* CVODE ODE solver parameters (default to type's quiet_NAN) */
  /* initial uniform thermodynamics */
  double P_INIT = std::numeric_limits<double>::signaling_NaN();    // initial pressure [Pa]
  double TEMP_INIT = std::numeric_limits<double>::signaling_NaN(); // initial parcel temperature [T]
  double relh_init = std::numeric_limits<double>::signaling_NaN(); // initial relative humidity (%)

  /* ODE parameters */
  double W_AVG = std::numeric_limits<double>::signaling_NaN();        // average amplitude of sinusoidal vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
  double T_HALF = std::numeric_limits<double>::signaling_NaN();       // timescale for w sinusoid, tau_half = T_HALF/pi [s]
  double cvode_rtol = std::numeric_limits<double>::signaling_NaN();   // relative tolerance for [P, T, qv, qc] ODEs integration
  double cvode_atol = std::numeric_limits<double>::signaling_NaN();   // absolute tolerances for [P, T, qv, qc] ODEs integration

  Config(const std::string_view config_filename)
  /* set input paramters as members of config
  class instance from txt configuration file */
  {
    std::cout << "\n--- configuration ---\n";
    
    loadconfiguration(config_filename);

    /* copy setup (config and constants files) to a txt file */
    const std::string filestr(config_filename);
    copyfiles2txt(setuptxt, {filestr, constants_filename});
    
    std::cout << "--- configuration: success ---\n";
  };
};

#endif // CONFIG_HPP