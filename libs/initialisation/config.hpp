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
  /* initialisation and output parameters */
  std::string initSDs_filename;      // binary filename for initialisation of superdrops
  std::string grid_filename;         // binary filename for gridbox boundaries
  std::string setuptxt;              // text filename to copy inital setup to as output
  std::filesystem::path zarrbasedir; // zarr store base directory
  int maxchunk;                      // maximum no. of elements in chunks of zarr store array

  /* SDM timestepping parameters */
  int cond_maxiters;  // maximum no. iterations of Newton Raphson Method
  double cond_rtol;   // relative tolerance for implicit euler integration
  double cond_atol;   //  abolute tolerance for implicit euler integration
  double CONDTSTEP;   // time between SD condensation events [s]
  double COLLTSTEP;   // time between SD collision events [s]
  double MOTIONTSTEP; // time between SD coordinate position updates [s]
  double COUPLTSTEP;  // time between SDM data output and thermodynamic coupling [s]
  double T_END;       // time span of integration [s]

  /* superdroplet parameters */
  int nSDsvec;       // initial no. elements in SDsInGBxs vector
  int SDnspace;      // number of spatial coordinates of superdroplets
  bool wetradiiinit; // set initial SD radii to equilibrium wet radius

  /* read in thermodynamics file parameters */
  std::string press_filename; // binary filename for pressure
  std::string temp_filename;  // binary filename for temperature
  std::string qvap_filename;  // binary filename for vapour mixing ratio
  std::string qcond_filename; // binary filename for liquid mixing ratio
  std::string wvel_filename;  // binary filename for vertical (z) velocity
  std::string uvel_filename;  // binary filename for horizontal x velocity
  std::string vvel_filename;  // binary filename for horizontal y velocity

  // /* CVODE ODE solver paramters */
  // double P_INIT;       // initial pressure [Pa]
  // double TEMP_INIT;    // initial parcel temperature [T]
  // double relh_init;    // initial relative humidity (%)
  // double qc_init;      // initial liquid water content []
  
  // bool doCouple;       // enable coupling from SDM to CVODE
  // bool doThermo;       // enable condensational growth of superdroplets
  // double W_AVG;        // average amplitude of sinusoidal vertical parcel speed [m/s] (dP/dt ~ w*dP/dz)
  // double T_HALF;       // timescale for w sinusoid, tau_half = T_HALF/pi [s]
  // double cvode_rtol;   // relative tolerance (tol) for integration
  // double cvode_atol_p; // absolute tolerances for thermodynamics ODEs [P, T, qv, qc]
  // double cvode_atol_temp;
  // double cvode_atol_qv;
  // double cvode_atol_qc;

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