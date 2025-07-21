/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: init_all_supers_from_binary.cpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct for reading in all super-droplets' initial conditions for CLEO SDM
 * (e.g. superdroplet attributes) from a binary file. InitAllSupersFromBinary instance
 * can be used by InitConds struct as SuperdropInitConds type.
 */

#include "initialise/init_all_supers_from_binary.hpp"

/* check all the vectors in the initdata struct all have sizes consistent with one another
and with maxnsupers. Include coords data in check if nspacedims > 0 */
void check_initdata_sizes(const InitSupersData &in, const size_t maxnsupers,
                          const size_t nspacedims) {
  std::vector<size_t> sizes({maxnsupers, in.sdgbxindexes.size(), in.xis.size(), in.radii.size(),
                             in.msols.size(), in.sdIds.size()});

  switch (nspacedims) {
    case 3:  // 3-D model
      sizes.push_back(in.coord2s.size());
      [[fallthrough]];
    case 2:  // 3-D or 2-D model
      sizes.push_back(in.coord1s.size());
      [[fallthrough]];
    case 1:  // 3-D, 2-D or 1-D model
      sizes.push_back(in.coord3s.size());
  }

  check_vectorsizes(sizes);
}

/* sets initial data for solutes as
a single SoluteProprties instance */
void InitAllSupersFromBinary::initdata_for_solutes(InitSupersData &initdata) const {
  initdata.solutes.at(0) = SoluteProperties{};
}

/* sets initial data for sdIds using an sdId's generator */
void InitAllSupersFromBinary::initdata_for_sdIds(InitSupersData &initdata) const {
  auto sdIdgen = Superdrop::IDType::Gen();

  for (size_t kk(0); kk < maxnsupers; ++kk) {
    initdata.sdIds.push_back(sdIdgen.next());
  }
}

/* sets initial data in initdata using data read
from a binary file called initsupers_filename */
void InitAllSupersFromBinary::initdata_from_binary(InitSupersData &initdata) const {
  std::ifstream file(open_binary(initsupers_filename));

  std::vector<VarMetadata> meta(metadata_from_binary(file));

  read_initdata_binary(initdata, file, meta);

  file.close();
}

/* copy data for vectors from binary file to initdata struct */
void InitAllSupersFromBinary::read_initdata_binary(InitSupersData &initdata, std::ifstream &file,
                                                   const std::vector<VarMetadata> &meta) const {
  initdata.sdgbxindexes = vector_from_binary<unsigned int>(file, meta.at(0));

  initdata.xis = vector_from_binary<uint64_t>(file, meta.at(1));

  initdata.radii = vector_from_binary<double>(file, meta.at(2));

  initdata.msols = vector_from_binary<double>(file, meta.at(3));

  initdata.coord3s = vector_from_binary<double>(file, meta.at(4));

  initdata.coord1s = vector_from_binary<double>(file, meta.at(5));

  initdata.coord2s = vector_from_binary<double>(file, meta.at(6));
}

/* data size returned is number of variables as
declared by the metadata for the first variable
in the initsupers file */
size_t InitAllSupersFromBinary::fetch_data_size() const {
  std::ifstream file(open_binary(initsupers_filename));

  VarMetadata meta(metadata_from_binary(file).at(0));

  file.close();

  return meta.nvar;
}
