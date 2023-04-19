// Author: Clara Bayley
// File: thermodynamicsfromfile.cpp
/* implementation of things specifically
to run uncoupled SDM where thermodynamics
are read from file */

#include "thermodynamicsfromfile.hpp"

std::vector<double> thermodynamicvar_from_binary(std::string_view filename)
{
  /* open file and read in the metatdata
  for all the variables in that file */
  std::ifstream file(open_binary(filename));
  std::vector<VarMetadata> meta(metadata_from_binary(file));
  
  /* read in the data for the 1st variable in the file */
  std::vector<double>
      thermovar(vector_from_binary<double>(file, meta.at(0)));

  return thermovar;
}

ThermodynamicsFromFile::ThermodynamicsFromFile(const Config &config,
                                               const int ngridboxes)
  : press(thermodynamicvar_from_binary(config.press_filename)),
  temp(thermodynamicvar_from_binary(config.temp_filename)),
  qvap(thermodynamicvar_from_binary(config.press_filename)),
  qcond(thermodynamicvar_from_binary(config.press_filename)),
  wvel(), uvel(), vvel()
{
  std::string windstr = set_windvelocities(config);
  std::cout << "\nFinished reading thermodynamics from binaries for:\n"
               "  pressure,\n  temperature,\n"
               "  water vapour mass mixing ratio,\n"
               "  liquid water mass mixing ratio,\n  "
            << windstr << '\n';
}

std::string ThermodynamicsFromFile::set_windvelocities(const Config &config)
{
  const int SDnspace(config.SDnspace);

  if (SDnspace == 0)
  {
    return "0-D model has no wind data";
  }

  else if (SDnspace <= 3) // means 1 <= SDnspace < 4
  {
    std::string info(std::to_string(SDnspace)+"-D model ");
    wvel = thermodynamicvar_from_binary(config.wvel_filename);
    
    if (SDnspace >= 2)
    {
      uvel = thermodynamicvar_from_binary(config.uvel_filename);

      if (SDnspace == 3)
      {
        vvel = thermodynamicvar_from_binary(config.vvel_filename);
        return info+"[w, u, v] wind velocity";
      }      
      return info+"[w, u] wind velocity";
    }
    return info+"vertical w wind velocity";
  }

  else // means SDnspace > 3
  {
    const std::string errmsg("SDnspace > 3 is invalid");
    throw std::invalid_argument(errmsg);
  }
}

void ThermodynamicsFromFile::run_thermostep(const int couplstep) const
{
  std::cout << "thermostep\n";
}