/*
 * ----- CLEO -----
 * File: gbxbounds_frombinary.cpp
 * Project: initialise
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 1st November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions for reading gridbox boundaries from
 * a binary file (used to then create a
 * map from gbxindexes to gridbox boundaries
 * for CLEO SDM e.g. a CartesianMaps)
 */


#include "./gbxbounds_frombinary.hpp"

GridboxBoundsFromBinary::
    GridboxBoundsFromBinary(const unsigned int nspacedims,
                            std::string_view gridfile)
/* read metadata and data in binary file called 'gridfile', then
return GridBoxBoundaries instance created from that data */
{
  /* open file and read in the metatdata
  for all the variables in gridfile */
  std::ifstream file(open_binary(gridfile));
  std::vector<VarMetadata> meta(metadata_from_binary(file));

  std::vector<size_t>
      ndims(vector_from_binary<size_t>(file, meta.at(0)));

  std::vector<unsigned int>
      gbxidxs(vector_from_binary<unsigned int>(file, meta.at(1))); 

  std::vector<double>
      gbxbounds(vector_from_binary<double>(file, meta.at(2))); 

  file.close();

  if (6 * gbxidxs.size() != gbxbounds.size() && gbxbounds.size() >=6)
  {
    std::string errormsg = "sizes of gbxidxs and gbxbounds vectors"
                           " read from gridfile not consistent";
    throw std::invalid_argument(errormsg);
  }

  is_gridbounds_SDnspace_compatible(SDnspace, gbxbounds, ndims);

  return GridBoxBoundaries{ndims, gbxidxs, gbxbounds};
}
