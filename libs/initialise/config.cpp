/*
 * ----- CLEO -----
 * File: config.cpp
 * Project: initialise
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functions involved in reading values from config files using
 * Config class
 */


#include "config.hpp"

void Config::loadconfiguration(const std::string_view config_filename)
/* read configuration file given by config_filename
  then calls configvariable to assign value to corresponding
  member of Config structure */
{
  std::ifstream file;
  open_file(config_filename, file);

  std::string line;
  std::cout << "----- Reading config file contents -----\n";
  while (getline(file, line))
  {
    line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

    if (line[0] == '*' || line[0] == '#' || line[0] == '/' || line.empty())
    {
      continue;
    }
    else
    {
      NameValue nv = getvariable_fromline(line);
      configvariable(nv.name, nv.value);
    }
  }

  file.close();
  std::cout << "---- Contents has been read, config file closed -----\n";
}

void Config::open_file(const std::string_view config_filename,
                       std::ifstream &file) const
/* opens file or throws error if opening fails */
{
  const std::string filestr(config_filename);
  std::cout << "opening config file: " << filestr << '\n';
  file.open(filestr);

  if (!file.is_open())
  {
    throw std::invalid_argument("Cannot open "+ filestr);
  }
}

Config::NameValue Config::getvariable_fromline(const std::string line) const
/* extracts name of variable and it's value from a string
of the form 'variable_name = value    #comment' */
{
  const std::string::size_type delimiter_pos = line.find("=");
  const std::string::size_type end_pos = line.find("#");
  const std::string name = line.substr(0, delimiter_pos);
  const std::string value = line.substr(delimiter_pos + 1,
                                        end_pos - delimiter_pos - 1);

  return {name, value};
}

bool Config::string2bool(const std::string value) const
/* if string == 'true', 'True' or '1', returns true
C++ boolean, else returns false */
{
  if (value == "true" || value == "True" || value == "1")
  {
    return true;
  }
  else
  {
    return false;
  }
}

void Config::configvariable(const std::string name, const std::string value)
/* setter function. assigns value of member of Config struct
  called 'name' by coverting strings containing
  value into actual value for that  members's type. If
  'name' member cannot be assigned, throw error */
{
  bool issuccess = false;

  /*** Initialisation and Output Data parameters ***/
  if (name == "constants_filename")
  {
    constants_filename = value;
    issuccess = true;
  }
  else if (name == "initsupers_filename")
  {
    initsupers_filename = value;
    issuccess = true;
  }
  else if (name == "grid_filename")
  {
    grid_filename = value;
    issuccess = true;
  }
  else if (name == "setuptxt")
  {
    setuptxt = value;
    issuccess = true;
  }
  else if (name == "zarrbasedir")
  {
    zarrbasedir = value;
    issuccess = true;
  }
  else if (name == "maxchunk")
  {
    maxchunk = stoi(value);
    issuccess = true;
  }

  /*** SDM Runtime parameters ***/
  /* domain setup */
  else if (name == "nspacedims")
  {
    nspacedims = stoi(value);
    issuccess = true;
  }
  else if (name == "ngbxs")
  {
    ngbxs = stoi(value);
    issuccess = true;
  }
  else if (name == "totnsupers")
  {
    totnsupers = stoi(value);
    issuccess = true;
  }

  /* timestepping */
  else if (name == "CONDTSTEP")
  {
    CONDTSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "COLLTSTEP")
  {
    COLLTSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "MOTIONTSTEP")
  {
    MOTIONTSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "COUPLTSTEP")
  {
    COUPLTSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "OBSTSTEP")
  {
    OBSTSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "T_END")
  {
    T_END = stod(value);
    issuccess = true;
  }

  /* microphysics */
  else if (name == "cond_iters")
  {
    cond_iters = stoul(value);
    issuccess = true;
  }
  else if (name == "cond_SUBTSTEP")
  {
    cond_SUBTSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "cond_rtol")
  {
    cond_rtol = stod(value);
    issuccess = true;
  }
  else if (name == "cond_atol")
  {
    cond_atol = stod(value);
    issuccess = true;
  }
  else if (name == "nfrags")
  {
    nfrags = stod(value);
    issuccess = true;
  }

  /* superdroplets */
  else if (name == "wetradiiinit")
  {
    wetradiiinit = string2bool(value);
    issuccess = true;
  }
  else if (name == "doAlterThermo")
  {
    doAlterThermo = string2bool(value);
    issuccess = true;
  }

  /*** Coupled Dynamics Solver Parameters ***/
  else if (name == "thermosolver")
  {
    thermosolver = value;
    issuccess = true;
  }
  if (thermosolver == "fromfile")
  {
    /* read in dynamics from file */
    configvariable_thermosolverfromfile(name, value);
    issuccess = true;
  }
  else if (thermosolver == "cvode")
  {
    /* CVODE ODE solver parameters */
    configvariable_thermosolvercvode(name, value);
    issuccess = true;
  }

  if (issuccess)
  {
    std::cout << "assigned " << name << " = " << value << '\n';
  }
  else
  {
    throw std::invalid_argument(name + " cannot be assigned with input value");
  }
}

void Config::configvariable_thermosolverfromfile(const std::string name,
                                                 const std::string value)
/* setter function for assigning 'value' to members of
Config struct called 'name' specifically for members
involved when thermosolver == 'fromfile'. Returns false only if
name is one of these specific members and also cannot be assigned */
{
  if (name == "press_filename")
  {
    press_filename = value;
  }
  else if (name == "temp_filename")
  {
    temp_filename = value;
  }
  else if (name == "qvap_filename")
  {
    qvap_filename = value;
  }
  else if (name == "qcond_filename")
  {
    qcond_filename = value;
  }
  else if (name == "wvel_filename")
  {
    wvel_filename = value;
  }
  else if (name == "uvel_filename")
  {
    uvel_filename = value;
  }
  else if (name == "vvel_filename")
  {
    vvel_filename = value;
  }
}

void Config::configvariable_thermosolvercvode(const std::string name,
                                              const std::string value)
/* setter function for assigning 'value' to members of
Config struct called 'name' specifically for members
involved when thermosolver == 'cvode' */
{
  /* initial (uniform) thermodynamic conditions */
  if (name == "P_INIT")
  {
    P_INIT = stod(value);
  }
  else if (name == "TEMP_INIT")
  {
    TEMP_INIT = stod(value);
  }
  else if (name == "relh_init")
  {
    relh_init = stod(value);
  }

  /* ODE parameters */
  else if (name == "W_AVG")
  {
    W_AVG = stod(value);
  }
  else if (name == "T_HALF")
  {
    T_HALF = stod(value);
  }
  else if (name == "cvode_rtol")
  {
    cvode_rtol = stod(value);
  }
  else if (name == "cvode_atol")
  {
    cvode_atol = stod(value);
  }
}