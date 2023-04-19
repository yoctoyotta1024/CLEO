// Author: Clara Bayley
// File: thermodynamicsfromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#ifndef THERMODYNAMICSFROMFILE_HPP
#define THERMODYNAMICSFROMFILE_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <istream>
#include <string>
#include <string_view>

#include "initialisation/config.hpp"
#include "initialisation/readbinary.hpp"

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

  std::string set_windvelocities(const Config &config);
  
  void check_thermodyanmics_vectorsizes(const int SDnspace,
                                        const size_t sz0);

public:
  ThermodynamicsFromFile(const Config &config, const int ngridboxes);

  void run_thermostep(const int couplstep) const;
};

#endif // THERMODYNAMICSFROMFILE_HPP 
