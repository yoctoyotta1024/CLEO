/*
 * ----- CLEO -----
 * File: initsupers_frombinary.cpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 31st October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for superdroplets' initial conditions
 * for CLEO SDM (e.g. superdroplet attributes)
 * by reading binary file. InitSupersFromBinary 
 * instance can be used by InitConds
 * struct as SuperdropInitConds type
 */


#include "./initsupers_frombinary.hpp"

InitSupersData InitSupersFromBinary::fetch_data() const
/* return InitSupersData created by reading a binary
file and creating a SoluteProperties struct.
Then check that the input data has the correct sizes. */
{
  InitSupersData initdata;

  init_solutes_data(initdata);
  initdata_from_binary(initdata); 

  check_initdata_sizes(initdata);

  return initdata;
}

void InitSupersFromBinary::
    init_solutes_data(InitSupersData &initdata) const
/* sets initial data for solutes as
a single SoluteProprties instance */
{
  initdata.solutes.at(0) = SoluteProperties{};
}

void InitSupersFromBinary::
    initdata_from_binary(InitSupersData &initdata) const
/* sets initial data in initdata using data read
from a binary file called initsupers_filename */
{
  std::ifstream file(open_binary(initsupers_filename));

  std::vector<VarMetadata> meta(metadata_from_binary(file));

  read_initdata_binary(initdata, file, meta);

  file.close();
};

void InitSupersFromBinary::
    check_initdata_sizes(const InitSupersData &in) const
/* check all the vectors in the initdata struct all
have sizes consistent with one another. Include
coords data in check if nspacedims != 0 */    
{
  std::vector<size_t> sizes({in.sdgbxindexes.size(),
                             in.xis.size(),
                             in.radii.size(),
                             in.msols.size()});

  if (nspacedims > 0)
  {
    sizes.push_back(in.coord3s.size());

    if (nspacedims > 1)
    {
      sizes.push_back(in.coord1s.size()); 

      if (nspacedims == 3)
      {
        sizes.push_back(in.coord2s.size());
      }
    }
  }

  check_vectorsizes(sizes);
}

void InitSupersFromBinary::
    read_initdata_binary(InitSupersData &initdata,
                         std::ifstream &file,
                         const std::vector<VarMetadata> &meta) const
/* copy data for vectors from binary file to initdata struct */
{
  initdata.sdgbxindexes = vector_from_binary<unsigned int>(file, meta.at(0));
  
  initdata.xis = vector_from_binary<unsigned long long>(file, meta.at(1));

  initdata.radii = vector_from_binary<double>(file, meta.at(2));

  inidata.msols = vector_from_binary<double>(file, meta.at(3));

  initdata.coord3s = vector_from_binary<double>(file, meta.at(4));

  initdata.coord1s = vector_from_binary<double>(file, meta.at(5));
  
  initdata.coord2s = vector_from_binary<double>(file, meta.at(6));
}

size_t InitSupersFromBinary::fetch_data_size() const
/* data size returned is number of variables as
declared by the metadata for the first variable 
in the initsupers file */
{
  std::ifstream file(open_binary(initsupers_filename));
  
  VarMetadata meta(metadata_from_binary(file).at(0)); 
  
  file.close();
   
  return meta.nvar; //TODO
}