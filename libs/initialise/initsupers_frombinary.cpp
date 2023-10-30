/*
 * ----- CLEO -----
 * File: initsupers_frombinary.cpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 30th October 2023
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
{
  InitSupersData initdata;

  init_solutes_data(initdata);
  initdata_from_binary(initdata); 

  return initdata;
}

void InitSupersFromBinary::
    init_solutes_data(InitSupersData &initdata) const
{
  initdata.solutes.at(0) = SoluteProperties{};
}

void InitSupersFromBinary::
    initdata_from_binary(InitSupersData &initdata) const
{
  std::ifstream file(open_binary(initsupers_filename));

  std::vector<VarMetadata> meta(metadata_from_binary(file));

  read_initdata_binary(initdata, file, meta);

  file.close();

  // check_vectorsizes({sd_gbxindex.size(), eps.size(),
  //                   radius.size(), m_sol.size()});

  // return InitSDsData{sd_gbxindex, eps, radius, m_sol,
  //                    coord3, coord1, coord2};
}

void InitSupersFromBinary::
    read_initdata_binary(InitSupersData &initdata,
                         std::ifstream &file,
                         const VarMetadata &varmeta) const
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
{
  return 0; //TODO
}