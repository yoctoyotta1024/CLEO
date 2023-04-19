// Author: Clara Bayley
// File: thermodynamicsfromfile.cpp
/* implementation of things specifically
to run uncoupled SDM where thermodynamics
are read from file */

#include "thermodynamicsfromfile.hpp"

std::vector<double> thermovar_frombinary(std::string_view filename)
{
  std::ifstream file(open_binary(filename));

  std::vector<VarMetadata> meta(metadata_from_binary(file));
  
  VarMetadata varmeta(meta.at(0));
  file.seekg(varmeta.b0, std::ios::beg);
  std::vector<double> thermovar(varmeta.nvar, 0);
  binary_into_buffer<double>(file, thermovar);

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