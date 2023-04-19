// Author: Clara Bayley
// File: read_initsuperdrops.cpp
/* Functionality to create InitSDsData
struct by reading data from a binary file */

#include "read_initsuperdrops.hpp"

InitSDsData get_initsuperdropsdata(std::string_view initSDsfile)
{
  std::ifstream file(open_binary(initSDsfile));

  std::vector<VarMetadata> meta(metadata_from_binary(file));

  std::vector<unsigned int>
      sd_gbxindex(vector_from_binary<unsigned int>(file, meta.at(0)));
  
  std::vector<size_t>
      eps(vector_from_binary<size_t>(file, meta.at(1)));

  std::vector<double>
      radius(vector_from_binary<double>(file, meta.at(2)));

  std::vector<double>
      m_sol(vector_from_binary<double>(file, meta.at(3)));

  std::vector<double>
      coord3(vector_from_binary<double>(file, meta.at(4)));

  std::vector<double>
      coord1(vector_from_binary<double>(file, meta.at(5)));
  
  std::vector<double>
      coord2(vector_from_binary<double>(file, meta.at(6)));

  file.close();

  check_vectorsizes({sd_gbxindex.size(), eps.size(),
                    radius.size(), m_sol.size()});

  return InitSDsData{sd_gbxindex, eps, radius, m_sol,
                     coord3, coord1, coord2};
};