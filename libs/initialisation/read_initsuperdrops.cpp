// Author: Clara Bayley
// File: read_initsuperdrops.cpp
/* Functionality to create InitSDsData
struct by reading data from a binary file */

#include "read_initsuperdrops.hpp"

InitSDsData get_initsuperdropsdata(std::string_view initSDsfile)
{
  std::ifstream file(open_binary(initSDsfile));

  std::vector<VarMetadata> meta (metadata_from_binary(file));

  VarMetadata var(meta.at(0));
  file.seekg(var.b0, std::ios::beg);
  std::vector<unsigned int> isd_gbxindex(var.nvar, 0);
  binary_into_buffer<unsigned int>(file, isd_gbxindex);
 
  var = meta.at(1);
  file.seekg(var.b0, std::ios::beg);
  std::vector<size_t> ieps(var.nvar, 0);
  binary_into_buffer<size_t>(file, ieps);
 
  var = meta.at(2);
  file.seekg(var.b0, std::ios::beg);
  std::vector<double> iradius(var.nvar, 0);
  binary_into_buffer<double>(file, iradius);

  var = meta.at(3);
  file.seekg(var.b0, std::ios::beg);
  std::vector<double> im_sol(var.nvar, 0);
  binary_into_buffer<double>(file, im_sol);

  var = meta.at(4);
  file.seekg(var.b0, std::ios::beg);
  std::vector<double> icoord3(var.nvar, 0);
  binary_into_buffer<double>(file, icoord3);

  file.close();
  
  std::vector<double> icoord1(0);
  std::vector<double> icoord2(0);

  return InitSDsData{isd_gbxindex, ieps, iradius, im_sol,
                      icoord3, icoord1, icoord2};
};