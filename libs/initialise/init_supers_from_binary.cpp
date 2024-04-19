/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: init_supers_from_binary.cpp
 * Project: initialise
 * Created Date: Monday 30th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 19th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct for reading in some super-droplets' initial conditions for CLEO SDM
 * (e.g. superdroplet attributes) from a binary file. InitAllSupersFromBinary instance
 * can be used by InitConds struct as SuperdropInitConds type.
 */

#include "initialise/init_supers_from_binary.hpp"

template <typename T>
inline std::vector<T> nan_vector(const size_t size) {
  const auto nanValue = std::numeric_limits<T>::signaling_NaN();
  return std::vector<T>(size, nanValue);
}

/* sets sdIds for un-initialised superdrops' using an sdId's generator */
std::vector<Superdrop::IDType> InitSupersFromBinary::sdIds_for_uninitialised_superdrops(
    const size_t size) const {
  auto sdIdgen = Superdrop::IDType::Gen();

  auto sdIds = std::vector<Superdrop::IDType>(
      size, sdIdgen.set(std::numeric_limits<unsigned int>::signaling_NaN()));

  return sdIds;
}

/* adds data for un-initialised (and out of bounds) superdrops into initdata so that initial
conditions exist for maxnsupers number of superdrops in total */
InitSupersData InitSupersFromBinary::add_uninitialised_superdrops_data(
    InitSupersData &initdata) const {
  const auto size = maxnsupers - initdata.sdgbxindexes.size();

  const auto sdgbxindexes = std::vector<unsigned int>(size, LIMITVALUES::uintmax);  // out of bounds
  const auto coord3s = nan_vector<double>(size);
  const auto coord1s = nan_vector<double>(size);
  const auto coord2s = nan_vector<double>(size);
  const auto radii = nan_vector<double>(size);
  const auto msols = nan_vector<double>(size);
  const auto xis = nan_vector<uint64_t>(size);
  const auto sdIds = sdIds_for_uninitialised_superdrops(size);

  const auto nandata = InitSupersData{
      initdata.solutes, sdgbxindexes, coord3s, coord1s, coord2s, radii, msols, xis, sdIds};

  return initdata + nandata;
}
