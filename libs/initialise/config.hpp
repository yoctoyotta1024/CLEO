/*
 * ----- CLEO -----
 * File: config.hpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
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

#include <string>
#include <string_view>
#include <fstream>
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
  /* Initialisation Files and Output Data parameters */
  std::string constants_filename;    // filename containing values of physical constants
  std::string initSDs_filename;      // binary filename for initialisation of SDs
  std::string grid_filename;         // binary filename for GBx boundaries
  std::string setuptxt;              // name of .txt output file to copy setup to
  std::filesystem::path zarrbasedir; // zarr store base directory
  int maxchunk;                      // maximum no. of elements in chunks of zarr store array

  /* SDM parameters */
  /* timestepping parameters */
  unsigned int cond_iters;     // suggested no. iterations of Newton Raphson Method
  double cond_SUBTSTEP;        // smallest timestep in cases where substepping occurs [s]
  double cond_rtol;            // relative tolerance for implicit euler integration
  double cond_atol;            // abolute tolerance for implicit euler integration
  double nfrags;               // expected number of fragments in (real) breakup even, (use -1 for nan value)
  double CONDTSTEP;            // time between SD condensation events [s]
  double COLLTSTEP;            // time between SD collision events [s]
  double MOTIONTSTEP;          // time between SDM motion [s]
  double COUPLTSTEP;           // time between thermodynamic couplings [s]
  double OBSTSTEP;             // time between SDM observations [s]
  double T_END;                // time span of integration [s]

  /* SDs parameters */
  int totnsupers;           // initial no. elements in SDs' vector (=total initial no. of SDs)
  int nspacedims;           // no. of spatial dimensions to model
  bool wetradiiinit;        // set initial SD radii to equilibrium wet radius
  bool doAlterThermo;       // enable condensation to alter the thermodynamic state
  std::string thermosolver; // type of thermodynamic solver to configure

  /* Read in Thermodynamics File parameters (default to empty) */
  std::string press_filename = ""; // binary filename for pressure
  std::string temp_filename = "";  // binary filename for temperature
  std::string qvap_filename = "";  // binary filename for vapour mixing ratio
  std::string qcond_filename = ""; // binary filename for liquid mixing ratio
  std::string wvel_filename = "";  // binary filename for vertical (z) velocity
  std::string uvel_filename = "";  // binary filename for horizontal x velocity
  std::string vvel_filename = "";  // binary filename for horizontal y velocity

  /* CVODE ODE solver parameters (default to type's quiet_NAN) */
  /* initial (uniform) thermodynamic conditions */
  double P_INIT = std::numeric_limits<double>::signaling_NaN();    // initial pressure [Pa]
  double TEMP_INIT = std::numeric_limits<double>::signaling_NaN(); // initial parcel temperature [T]
  double relh_init = std::numeric_limits<double>::signaling_NaN(); // initial relative humidity (%)
  double qc_init = std::numeric_limits<double>::signaling_NaN();   // initial liquid water content []

  /* ODE parameters */
  bool doThermo = std::numeric_limits<bool>::signaling_NaN();         // enable ODEs for adiabatic expansion
  double W_AVG = std::numeric_limits<double>::signaling_NaN();        // average amplitude of sinusoidal vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
  double T_HALF = std::numeric_limits<double>::signaling_NaN();       // timescale for w sinusoid, tau_half = T_HALF/pi [s]
  double cvode_rtol = std::numeric_limits<double>::signaling_NaN();   // relative tolerance for [P, T, qv, qc] ODEs integration
  double cvode_atol_p = std::numeric_limits<double>::signaling_NaN(); // absolute tolerances for [P, T, qv, qc] ODEs integration
  double cvode_atol_temp = std::numeric_limits<double>::signaling_NaN();
  double cvode_atol_qv = std::numeric_limits<double>::signaling_NaN();
  double cvode_atol_qc = std::numeric_limits<double>::signaling_NaN();

  Config(const std::string_view config_filename)
  /* set input paramters as members of config
  class instance from txt configuration file */
  {
    loadconfiguration(config_filename);

    /* copy setup (config and constants files) to a txt file */
    const std::string filestr(config_filename);
    CopyFiles2Txt()(setuptxt, {filestr, constants_filename});
  };
};

#endif // CONFIG_HPP