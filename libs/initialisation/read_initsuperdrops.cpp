// Author: Clara Bayley
// File: read_initsuperdrops.cpp
/* Functionality to create InitSDsData
struct by reading data from a binary file */

#include "read_initsuperdrops.hpp"

void check_vectorsizes(const std::vector<size_t> &sizes);
/* raise error if values in sizes vector are not the same. Used to check
if vectors for SD attributes created from reading initSDsfile and used to
male InitSdsData object are the same size */

InitSDsData get_initsuperdropsdata(std::string_view initSDsfile)
{
  std::ifstream file(open_binary(initSDsfile));

  std::vector<VarMetadata> meta(metadata_from_binary(file));

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

  var = meta.at(5);
  file.seekg(var.b0, std::ios::beg);
  std::vector<double> icoord1(var.nvar, 0);
  binary_into_buffer<double>(file, icoord1);

  var = meta.at(6);
  file.seekg(var.b0, std::ios::beg);
  std::vector<double> icoord2(var.nvar, 0);
  binary_into_buffer<double>(file, icoord2);

  file.close();

  check_vectorsizes({isd_gbxindex.size(), ieps.size(),
                   iradius.size(), im_sol.size()});

  return InitSDsData{isd_gbxindex, ieps, iradius, im_sol,
                     icoord3, icoord1, icoord2};
};

void check_vectorsizes(const std::vector<size_t> &sizes)
/* raise error if values in sizes vector are not the same. Used to check
if vectors for SD attributes created from reading initSDsfile and used to
male InitSdsData object are the same size */
{
  const size_t sz0 = sizes.front();
  for (auto sz : sizes)
  {
    if (sz != sz0)
    {
      const std::string err("sizes of vectors for InitSDsData"
                            "are not consistent");
      throw std::invalid_argument(err);
    }
  }
}