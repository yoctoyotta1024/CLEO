// Author: Clara Bayley
// File: thermodynamicsfromfile.cpp
/* implementation of things specifically
to run uncoupled SDM where thermodynamics
are read from file */

#include "thermodynamicsfromfile.hpp"

std::vector<double> thermovar_frombinary(std::string_view filename)
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
{
  thermovar_frombinary(config.press_filename);
  std::cout << "Pressure data read from binary\n";

}

void ThermodynamicsFromFile::run_thermostep(const int couplstep) const
{
  std::cout << "thermostep\n";
}