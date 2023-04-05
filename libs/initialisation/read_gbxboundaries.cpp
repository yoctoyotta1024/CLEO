// Author: Clara Bayley
// File: read_gbxboundaries.cpp
/* initialisatin of Maps4GridBoxes 
struct from binary file */

#include "read_gbxboundaries.hpp"

GridBoxBoundaries read_gbxboundaries(std::string_view gridfile)
/* read metadata and data in binary file called 'gridfile', then
return GridBoxBoundaries instance created from that data */
{
  std::ifstream file(open_binary(gridfile));

  std::vector<VarMetadata> meta(metadata_from_binary(file));
  
  VarMetadata var(meta.at(0));
  file.seekg(var.b0, std::ios::beg);
  std::vector<unsigned int> gbxidxs(var.nvar, 0);
  binary_into_buffer<unsigned int>(file, gbxidxs);
 
  var = meta.at(1);
  file.seekg(var.b0, std::ios::beg);
  std::vector<double> gbxbounds(var.nvar, 0);
  binary_into_buffer<double>(file, gbxbounds);

  file.close();

  return GridBoxBoundaries{gbxidxs, gbxbounds};
}