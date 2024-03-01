/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: config.cpp
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
 * File Description:
 * Functions involved in reading values from config files using
 * Config class
 */

#include "initialise/config.hpp"

/* read configuration file given by config_filename
  then calls configvariable to assign value to corresponding
  member of Config structure */
void Config::loadconfiguration(const std::string_view config_filename) {
  std::ifstream file;
  open_file(config_filename, file);

  std::string line;
  std::cout << "----- Reading config file contents -----\n";
  while (getline(file, line)) {
    line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

    if (line[0] == '*' || line[0] == '#' || line[0] == '/' || line.empty()) {
      continue;
    } else {
      NameValue nv = getvariable_fromline(line);
      configvariable(nv.name, nv.value);
    }
  }

  file.close();
  std::cout << "---- Contents has been read, config file closed -----\n";
}

/* opens file or throws error if opening fails */
void Config::open_file(const std::string_view config_filename, std::ifstream &file) const {
  const std::string filestr(config_filename);
  std::cout << "opening config file: " << filestr << '\n';
  file.open(filestr);

  if (!file.is_open()) {
    throw std::invalid_argument("Cannot open " + filestr);
  }
}

/* extracts name of variable and it's value from a string
of the form 'variable_name = value    #comment' */
Config::NameValue Config::getvariable_fromline(const std::string line) const {
  const std::string::size_type delimiter_pos = line.find("=");
  const std::string::size_type end_pos = line.find("#");
  const std::string name = line.substr(0, delimiter_pos);
  const std::string value = line.substr(delimiter_pos + 1, end_pos - delimiter_pos - 1);

  return {name, value};
}

/* if string == 'true', 'True' or '1', returns true
C++ boolean, else returns false */
bool Config::string2bool(const std::string value) const {
  if (value == "true" || value == "True" || value == "1") {
    return true;
  } else {
    return false;
  }
}

/* setter function. assigns value of member of Config struct
  called 'name' by coverting strings containing
  value into actual value for that  members's type. If
  'name' member cannot be assigned, throw error */
void Config::configvariable(const std::string name, const std::string value) {
  bool issuccess = false;

  /*** Initialisation and Output Data parameters ***/
  if (name == "constants_filename") {
    constants_filename = value;
    issuccess = true;
  } else if (name == "initsupers_filename") {
    initsupers_filename = value;
    issuccess = true;
  } else if (name == "grid_filename") {
    grid_filename = value;
    issuccess = true;
  } else if (name == "setuptxt") {
    setuptxt = value;
    issuccess = true;
  } else if (name == "stats_filename") {
    stats_filename = value;
    issuccess = true;
  } else if (name == "zarrbasedir") {
    zarrbasedir = value;
    issuccess = true;
  } else if (name == "maxchunk") {
    maxchunk = stoi(value);
    issuccess = true;
  } else if (name == "nspacedims") {  /*** SDM Runtime parameters ***/
    nspacedims = stoi(value);         /* domain setup */
    issuccess = true;
  } else if (name == "ngbxs") {
    ngbxs = stoi(value);
    issuccess = true;
  } else if (name == "totnsupers") {
    totnsupers = stoi(value);
    issuccess = true;
  } else if (name == "CONDTSTEP") {   /* timestepping */
    CONDTSTEP = stod(value);
    issuccess = true;
  } else if (name == "COLLTSTEP") {
    COLLTSTEP = stod(value);
    issuccess = true;
  } else if (name == "MOTIONTSTEP") {
    MOTIONTSTEP = stod(value);
    issuccess = true;
  } else if (name == "COUPLTSTEP") {
    COUPLTSTEP = stod(value);
    issuccess = true;
  } else if (name == "OBSTSTEP") {
    OBSTSTEP = stod(value);
    issuccess = true;
  } else if (name == "T_END") {
    T_END = stod(value);
    issuccess = true;
  } else if (name == "cond_iters") {  /* microphysics */
    cond_iters = stoul(value);
    issuccess = true;
  } else if (name == "cond_SUBTSTEP") {
    cond_SUBTSTEP = stod(value);
    issuccess = true;
  } else if (name == "cond_rtol") {
    cond_rtol = stod(value);
    issuccess = true;
  } else if (name == "cond_atol") {
    cond_atol = stod(value);
    issuccess = true;
  } else if (name == "doAlterThermo") {  /* superdroplets */
    doAlterThermo = string2bool(value);
    issuccess = true;
  } else if (name == "thermosolver") {  /*** Coupled Dynamics Solver Parameters ***/
    thermosolver = value;
    issuccess = true;
  }

  if (thermosolver == "fromfile") {
    /* read in dynamics from file */
    configvariable_thermosolverfromfile(name, value);
    issuccess = true;
  } else if (thermosolver == "cvode") {
    /* CVODE ODE solver parameters */
    configvariable_thermosolvercvode(name, value);
    issuccess = true;
  }

  if (issuccess) {
    std::cout << "assigned " << name << " = " << value << '\n';
  } else {
    throw std::invalid_argument(name + " cannot be assigned with input value");
  }
}

/* setter function for assigning 'value' to members of
Config struct called 'name' specifically for members
involved when thermosolver == 'fromfile'. Returns false only if
name is one of these specific members and also cannot be assigned */
void Config::configvariable_thermosolverfromfile(const std::string name, const std::string value) {
  if (name == "press_filename") {
    press_filename = value;
  } else if (name == "temp_filename") {
    temp_filename = value;
  } else if (name == "qvap_filename") {
    qvap_filename = value;
  } else if (name == "qcond_filename") {
    qcond_filename = value;
  } else if (name == "wvel_filename") {
    wvel_filename = value;
  } else if (name == "uvel_filename") {
    uvel_filename = value;
  } else if (name == "vvel_filename") {
    vvel_filename = value;
  }
}

/* setter function for assigning 'value' to members of
Config struct called 'name' specifically for members
involved when thermosolver == 'cvode' */
void Config::configvariable_thermosolvercvode(const std::string name, const std::string value) {
  /* initial (uniform) thermodynamic conditions */
  if (name == "P_INIT") {
    P_INIT = stod(value);
  } else if (name == "TEMP_INIT") {
    TEMP_INIT = stod(value);
  } else if (name == "relh_init") {
    relh_init = stod(value);
  } else if (name == "W_AVG") {  /* ODE parameters */
    W_AVG = stod(value);
  } else if (name == "T_HALF") {
    T_HALF = stod(value);
  } else if (name == "cvode_rtol") {
    cvode_rtol = stod(value);
  } else if (name == "cvode_atol") {
    cvode_atol = stod(value);
  }
}
