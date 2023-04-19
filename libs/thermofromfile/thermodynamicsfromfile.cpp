// Author: Clara Bayley
// File: thermodynamicsfromfile.cpp
/* implementation of things specifically
to run uncoupled SDM where thermodynamics
are read from file */

#include "thermodynamicsfromfile.hpp"

ThermodynamicsFromFile::ThermodynamicsFromFile(const Config &config,
                                               const int ngridboxes)
{
  std::cout << "initit\n";
}

void ThermodynamicsFromFile::run_thermostep(const int couplstep) const
{
  std::cout << "thermostep\n";
}