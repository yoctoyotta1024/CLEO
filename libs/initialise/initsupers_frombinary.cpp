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

  init_solutes(initdata);
  initdata_from_binary(initdata); 

  return initdata;
}

void InitSupersFromBinary::
    init_solutes(InitSupersData &initdata) const
{
  initdata.solutes.at(0) = SoluteProperties{};
}

void InitSupersFromBinary::
    initdata_from_binary(InitSupersData &initdata) const
{
  // std::ifstream file(open_binary(initsupers_filename));

  // std::vector<VarMetadata> meta(metadata_from_binary(file));

  // std::vector<unsigned int>
  //     sd_gbxindex(vector_from_binary<unsigned int>(file, meta.at(0)));
  
  // std::vector<unsigned long long>
  //     eps(vector_from_binary<unsigned long long>(file, meta.at(1)));

  // std::vector<double>
  //     radius(vector_from_binary<double>(file, meta.at(2)));

  // std::vector<double>
  //     m_sol(vector_from_binary<double>(file, meta.at(3)));

  // std::vector<double>
  //     coord3(vector_from_binary<double>(file, meta.at(4)));

  // std::vector<double>
  //     coord1(vector_from_binary<double>(file, meta.at(5)));
  
  // std::vector<double>
  //     coord2(vector_from_binary<double>(file, meta.at(6)));

  // file.close();

  // check_vectorsizes({sd_gbxindex.size(), eps.size(),
  //                   radius.size(), m_sol.size()});

  // return InitSDsData{sd_gbxindex, eps, radius, m_sol,
  //                    coord3, coord1, coord2};
};

size_t InitSupersFromBinary::fetch_data_size() const
{
  return 0; //TODO
}