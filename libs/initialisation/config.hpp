// Author: Clara Bayley
// File: config.hpp
/* Header file for configuration Class
including functions involved
in reading values from config files */

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <filesystem>

#include "./copyfiles2txt.hpp"

struct Config
{
private:
  struct NameValue
  {
    std::string name;
    std::string value;
  };

  void loadconfiguration(const std::string configfilepath);
  /* read configuration file given by configfilepath
  then calls configvariable to assign value to corresponding
  member of Config structure */

  void open_file(const std::string configfilepath,
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
  value into actual value for that  members's type. If
  'name' member cannot be assigned, throw error */

public:
  /* in/output data parameters */
  std::string initSDs_filename; // binary filename for initialisation of superdrops
  std::string grid_filename;    // binary filename for gridbox boundaries
  std::string setuptxt;         // text filename to copy inital setup to as output
  std::filesystem::path zarrbasedir; // zarr store base directory
  int maxcsize;              // size of array chunks for SD attributes data if outputing to zarr store

  /* Model Settings */
  /* Model timestep parameters */
  int cond_maxiters;  // maximum no. iterations of Newton Raphson Method
  double cond_rtol;   // relative tolerance for implicit euler integration
  double cond_atol;   //  abolute tolerance for implicit euler integration
  double COND_TSTEP;  // time between SD condensation events [s]
  double COLL_TSTEP;  // time between SD collision events [s]
  double SEDI_TSTEP; // time between SD coordinate position updates [s]
  double XCHANGE_TSTEP; // time between SD exchange between gridboxes [s]
  double OUT_TSTEP;   // time between outputting data (and coupling) [s]
  double TEND;        // time span of integration [s]

  /* Superdroplet init params */
  int SDnspace;    // number of spatial coordinates of superdroplets
  int NSUPERS;     // max. no. distinct superdrop objects in array

  /* initial parcel conditions */
  double P_INIT;    // initial pressure [Pa]
  double TEMP_INIT; // initial parcel temperature [T]
  double relh_init; // initial relative humidity (%)
  double qc_init;   // initial liquid water content []

  /* CVODE ODE solver paramters */
  bool doCouple;       // enable coupling from SDM to CVODE
  bool doThermo;       // enable condensational growth of superdroplets
  double W_AVG;        // average amplitude of sinusoidal vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
  double T_HALF;       // timescale for w sinusoid, tau_half = T_HALF/pi [s]
  double cvode_rtol;   // relative tolerance (tol) for integration
  double cvode_atol_p; // absolute tolerances for thermodynamics ODEs [P, T, qv, qc]
  double cvode_atol_temp;
  double cvode_atol_qv;
  double cvode_atol_qc;

  Config(const std::string configfilepath)
  /* set input paramters as members of config
  class instance from txt configuration file */
  {
    loadconfiguration(configfilepath);
  };

  Config(const std::string configfilepath,
         const std::string constantsfilepath)
  /* set input paramters as members of config class instance
  from txt configuration file. Then also copy contents of
  configuration and constants files into 'setuptxt' file */
  {
    loadconfiguration(configfilepath);

    /* copy setup (config and constants files) to a txt file */
    CopyFiles2Txt()(setuptxt, {configfilepath, constantsfilepath});
  };
};

#endif // CONFIG_HPP