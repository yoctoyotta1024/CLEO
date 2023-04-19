// Author: Clara Bayley
// File: thermodynamicsfromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#ifndef THERMODYNAMICSFROMFILE_HPP
#define THERMODYNAMICSFROMFILE_HPP

#include <iostream>
#include <vector>

#include "initialisation/config.hpp"

class ThermodynamicsFromFile
{
private:
  std::vector<double> press;
  std::vector<double> temp;
  std::vector<double> qvap;
  std::vector<double> qcond;
  std::vector<double> wvel;
  std::vector<double> uvel;
  std::vector<double> vvel;

public:
  ThermodynamicsFromFile(const Config &config, const int ngridboxes);

  void run_thermostep(const int couplstep) const;
};

#endif // THERMODYNAMICSFROMFILE_HPP 
