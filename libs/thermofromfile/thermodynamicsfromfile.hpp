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
#include <functional>
#include <array>
#include <utility>

#include "initialisation/config.hpp"
#include "initialisation/readbinary.hpp"

class ThermodynamicsFromFile
{
private:
  const std::array<size_t, 3> ndims; // number of (centres of) gridboxes in [z,x,y] directions
  size_t atpos;                // position of thermodata for 0th gridbox at given timestep
  std::vector<double> press;
  std::vector<double> temp;
  std::vector<double> qvap;
  std::vector<double> qcond;
  std::vector<double> wvel; // w velocity define of z faces of gridboxes
  std::vector<double> uvel; // u velocity define of x faces of gridboxes
  std::vector<double> vvel; // v velocity define of y faces of gridboxes

  std::string set_windvelocities(const Config &config);

  void check_thermodyanmics_vectorsizes(const int SDnspace,
                                     const std::array<size_t, 3> &ndims,
                                     const size_t nsteps) const;

  double get_var(const std::vector<double> &vec,
                 const size_t ii) const
  /* use atpos to get a value from 'vec' (which should
  correspond to the thermodynamic variable for
  ii'th gridbox at desired timestep (step = atos//ngridboxes) */
  {
    return vec.at(atpos+ii);
  }

public:
  std::function<std::pair<double,double>(const unsigned int)> get_wvelzfaces; // funcs to get velocity defined in construction of class 
  std::function<std::pair<double,double>(const unsigned int)> get_uvelxfaces; // warning: these functions are not const member funcs by default
  std::function<std::pair<double,double>(const unsigned int)> get_vvelyfaces;

  ThermodynamicsFromFile(const Config &config,
                         const std::array<size_t, 3> &ndims,
                         const size_t nsteps);

  void run_thermostep()
  /* increment position of thermodata for 0th gridbox to positon
  at next timestep (ie. ngridboxes further along vector) */
  {
    atpos += (ndims[0]*ndims[1]*ndims[2]);
  }

  double get_press(const unsigned int gbxindex) const
  {
    return get_var(press, (size_t)gbxindex);
  }

  double get_temp(const unsigned int gbxindex) const
  {
    return get_var(temp, (size_t)gbxindex);
  }
  
  double get_qvap(const unsigned int gbxindex) const
  {
    return get_var(qvap, (size_t)gbxindex);
  }

  double get_qcond(const unsigned int gbxindex) const
  {
    return get_var(qcond, (size_t)gbxindex);
  }
};

#endif // THERMODYNAMICSFROMFILE_HPP 
