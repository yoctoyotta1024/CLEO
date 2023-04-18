// Author: Clara Bayley
// File: thermodynamicsfromfile.hpp
/* Header file for functions specifically
to run uncoupled SDM where thermodynamics
are read from file */

#ifndef THERMODYNAMICSFROMFILE_HPP
#define THERMODYNAMICSFROMFILE_HPP

//#include "initialisation/config.hpp"
#include <iostream>

class ThermodynamicsFromFile
{
  //ThermodynamicsFromFile(const Config &config){}
  ThermodynamicsFromFile()
  {
    std::cout << "here thermo from file is init-ed\n";
  }

  void run_thermostep(const int couplstep)
  {
    std::cout << "thermostep\n";
  }
}

#endif // THERMODYNAMICSFROMFILE_HPP 
