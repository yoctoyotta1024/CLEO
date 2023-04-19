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
  size_t atpos; // position of thermodata for 0th gridbox at given timestep
  size_t ngrid; // number of gridboxes in domain (=increment of atpos)  
  std::vector<double> press;
  std::vector<double> temp;
  std::vector<double> qvap;
  std::vector<double> qcond;
  std::vector<double> wvel;
  std::vector<double> uvel;
  std::vector<double> vvel;

  std::string set_windvelocities(const Config &config);
  
  void check_thermodyanmics_vectorsizes(const int SDnspace,
                                        const size_t sz0) const;

public:
  ThermodynamicsFromFile(const Config &config,
                         const size_t nsteps, const size_t ngridboxes);

  void run_thermostep()
  /* increment position of thermodata for 0th gridbox to positon
  at next timestep (ie. ngridboxes further along vector) */
  {
    atpos += ngrid;
  }

  double get_press(const unsigned int gbxindex) const
  {
    return press.at(atpos+(size_t)gbxindex);
  }

  double get_temp(const unsigned int gbxindex) const
  {
    return temp.at(atpos+(size_t)gbxindex);
  }
  
  double get_qvap(const unsigned int gbxindex) const
  {
    return qvap.at(atpos+(size_t)gbxindex);
  }

  double get_qcond(const unsigned int gbxindex) const
  {
    return qcond.at(atpos+(size_t)gbxindex);
  }

  double get_wvel(const unsigned int gbxindex) const
  {
    return wvel.at(atpos+(size_t)gbxindex);
  }

  double get_uvel(const unsigned int gbxindex) const
  {
    return uvel.at(atpos+(size_t)gbxindex);
  }

  double get_vvel(const unsigned int gbxindex) const
  {
    return vvel.at(atpos+(size_t)gbxindex);
  }
};

#endif // THERMODYNAMICSFROMFILE_HPP 
