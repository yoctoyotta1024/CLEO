// Author: Clara Bayley
// File: read_initsuperdrops.hpp
/* Header for initialisation of superdrops via
creation of an InitSDsData struct with data 
from a binary file */

#ifndef READ_INITSUPERDROPS_HPP
#define READ_INITSUPERDROPS_HPP

#include <vector>
#include <string>
#include <string_view>
#include <fstream>
#include <iostream>
#include <utility>

#include "readbinary.hpp"

struct InitSDsData
{
  std::vector<size_t> eps_init;
  std::vector<double> radius_init;
  std::vector<double> m_sol_init;
  std::vector<double> coord3_init; 
  std::vector<double> coord1_init;
  std::vector<double> coord2_init;
};

InitSDsData get_initsuperdropsdata(std::string_view initSDs_filename);

#endif // READ_INITSUPERDROPS_HPP 