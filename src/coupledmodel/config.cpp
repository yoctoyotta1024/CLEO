// Author: Clara Bayley
// File: config.cpp
/* Functions involved in reading
values from config files using
Config class */

#include "config.hpp"

void Config::loadconfiguration(const std::string configfilepath)
/* read configuration file given by configfilepath
  then calls configvariable to assign value to corresponding
  member of Config structure */
{
  std::ifstream file;
  open_file(configfilepath, file);

  std::string line;
  std::cout << "----- Reading config file contents -----\n";
  while (getline(file, line))
  {
    line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

    if (line[0] == '#' || line[0] == '/' || line.empty())
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

void Config::open_file(const std::string configfilepath,
                       std::ifstream &file) const
/* opens file or throws error if opening fails */
{
  std::cout << "opening config file: " << configfilepath << '\n';
  file.open(configfilepath);

  if (!file.is_open())
  {
    throw std::invalid_argument("Cannot open " + configfilepath);
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

void Config::configvariable(const std::string name, std::string value)
/* setter function. assigns value of member of Config struct
  called 'name' by coverting strings containing
  value into actual value for that  members's type. If
  'name' member cannot be assigned, throw error */
{
  bool issuccess = false;

  if (name == "initSDs_filename")
  {
    initSDs_filename = value;
    issuccess = true;
  }
  else if (name == "grid_filename")
  {
    grid_filename = value;
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
  else if (name == "maxcsize")
  {
    maxcsize = stoi(value);
    issuccess = true;
  }
  else if (name == "cond_maxiters")
  {
    cond_maxiters = stoi(value);
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
  else if (name == "COND_TSTEP")
  {
    COND_TSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "COLL_TSTEP")
  {
    COLL_TSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "SEDI_TSTEP")
  {
    SEDI_TSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "XCHANGE_TSTEP")
  {
    XCHANGE_TSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "OUT_TSTEP")
  {
    OUT_TSTEP = stod(value);
    issuccess = true;
  }
  else if (name == "TEND")
  {
    TEND = stod(value);
    issuccess = true;
  }

  else if (name == "SDnspace")
  {
    SDnspace = stoi(value);
    issuccess = true;
  }
  else if (name == "NSUPERS")
  {
    NSUPERS = stoi(value);
    issuccess = true;
  }

  else if (name == "P_INIT")
  {
    P_INIT = stod(value);
    issuccess = true;
  }
  else if (name == "TEMP_INIT")
  {
    TEMP_INIT = stod(value);
    issuccess = true;
  }
  else if (name == "relh_init")
  {
    relh_init = stod(value);
    issuccess = true;
  }
  else if (name == "qc_init")
  {
    qc_init = stod(value);
    issuccess = true;
  }

  else if (name == "doCouple")
  {
    doCouple = string2bool(value);
    issuccess = true;
  }
  else if (name == "doThermo")
  {
    doThermo = string2bool(value);
    issuccess = true;
  }
  else if (name == "W_AVG")
  {
    W_AVG = stod(value);
    issuccess = true;
  }
  else if (name == "T_HALF")
  {
    T_HALF = stod(value);
    issuccess = true;
  }
  else if (name == "cvode_rtol")
  {
    cvode_rtol = stod(value);
    issuccess = true;
  }
  else if (name == "cvode_atol_p")
  {
    cvode_atol_p = stod(value);
    issuccess = true;
  }
  else if (name == "cvode_atol_temp")
  {
    cvode_atol_temp = stod(value);
    issuccess = true;
  }
  else if (name == "cvode_atol_qv")
  {
    cvode_atol_qv = stod(value);
    issuccess = true;
  }
  else if (name == "cvode_atol_qc")
  {
    cvode_atol_qc = stod(value);
    issuccess = true;
  }

  if (issuccess)
  {
    std::cout << "assigned variable: " << name << " = " << value << '\n';
  }
  else
  {
    throw std::invalid_argument(name + " cannot be assigned with input value");
  }
}