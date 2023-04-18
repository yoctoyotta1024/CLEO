// Author: Clara Bayley
// File: thermodynamicsfromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#ifndef THERMODYNAMICSFROMFILE_HPP
#define THERMODYNAMICSFROMFILE_HPP

//#include "initialisation/config.hpp"
#include <iostream>
#include <vector>

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
  //ThermodynamicsFromFile(const Config &config){}
  ThermodynamicsFromFile(const int ngridboxes)
  {
    std::cout << "here thermo from file is init-ed eg. ngridboxes = "
    << ngridboxes << '\n';
  }

  void run_thermostep(const int couplstep) const;
};

#endif // THERMODYNAMICSFROMFILE_HPP 
