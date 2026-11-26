/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: init_supers_from_binary.hpp
 * Project: initialise
 * Created Date: Friday 19th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct for reading in some super-droplets' initial conditions for CLEO SDM
 * (e.g. superdroplet attributes) from a binary file. InitAllSupersFromBinary instance
 * can be used by InitConds struct as SuperdropInitConds type.
 */

#ifndef LIBS_INITIALISE_INIT_SUPERS_FROM_BINARY_HPP_
#define LIBS_INITIALISE_INIT_SUPERS_FROM_BINARY_HPP_

#include <mpi.h>

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "../cleoconstants.hpp"
#include "configuration/communicator.hpp"
#include "configuration/optional_config_params.hpp"
#include "initialise/init_all_supers_from_binary.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/readbinary.hpp"
#include "superdrops/superdrop.hpp"

template <typename T>
inline std::vector<T> init_nan_vector(const size_t size) {
  const auto nanValue = std::numeric_limits<T>::signaling_NaN();
  return std::vector<T>(size, nanValue);
}

/* struct containing functions which return data for the initial conditions needed to create
superdroplets e.g. via the CreateSupers struct */
template <typename GbxMaps>
struct InitSupersFromBinary {
 private:
  size_t maxnsupers;  /**< total number of super-droplets (in kokkos view on device) */
  size_t initnsupers; /**< initial no. of super-droplets to initialise */
  std::filesystem::path initsupers_filename; /**< filename for super-droplets' initial conditons */
  unsigned int nspacedims; /**< number of spatial dimensions to model (0-D, 1-D, 2-D of 3-D) */
  const GbxMaps& gbxmaps;  /**< hook to get to gridbox maps for current cartesian domain */

  /* returns InitSupersData created by reading some data from a binary file and
  filling the rest with un-initialised super-droplets */
  InitSupersData fetch_superdrops_from_file() const {
    auto initsupers = InitAllSupersFromBinary(initnsupers, initsupers_filename, nspacedims);
    return initsupers.fetch_data();
  }

  void trim_nonlocal_superdrops(InitSupersData& initdata) const {
    int my_rank;
    my_rank = init_communicator::get_comm_rank();

    if (gbxmaps.get_total_global_ngridboxes() == gbxmaps.get_local_ngridboxes_hostcopy()) return;

    // Go through all superdrops and resets the values of the non-local ones
    auto gbxindex = size_t{0};
    for (size_t superdrop_index = 0; superdrop_index < initdata.sdgbxindexes.size();) {
      gbxindex = initdata.sdgbxindexes[superdrop_index];
      if (my_rank != gbxmaps.get_domain_decomposition().get_gridbox_owner_process(gbxindex)) {
        // resets superdrops which are in gridboxes not owned by this process
        initdata.sdgbxindexes[superdrop_index] = LIMITVALUES::oob_gbxindex;
        initdata.xis[superdrop_index] = std::numeric_limits<uint64_t>::signaling_NaN();

        initdata.radii[superdrop_index] = std::numeric_limits<double>::signaling_NaN();

        initdata.msols[superdrop_index] = std::numeric_limits<double>::signaling_NaN();

        initdata.coord3s[superdrop_index] = std::numeric_limits<double>::signaling_NaN();

        initdata.coord1s[superdrop_index] = std::numeric_limits<double>::signaling_NaN();

        initdata.coord2s[superdrop_index] = std::numeric_limits<double>::signaling_NaN();

      } else {
        // updates superdrop gridbox index from global to local
        initdata.sdgbxindexes[superdrop_index] = gbxmaps.global_to_local_gbxindex(gbxindex);
      }
      superdrop_index++;
    }
  }

  /* adds data for un-initialised (and out of bounds) superdrops into initdata so that initial
  conditions exist for maxnsupers number of superdrops in total */
  InitSupersData add_uninitialised_superdrops_data(InitSupersData& initdata) const {
    const auto size = maxnsupers - initdata.sdgbxindexes.size();

    const auto sdgbxindexes =
        std::vector<unsigned int>(size, LIMITVALUES::oob_gbxindex);  // out of bounds
    const auto coord3s = init_nan_vector<double>(size);
    const auto coord1s = init_nan_vector<double>(size);
    const auto coord2s = init_nan_vector<double>(size);
    const auto radii = init_nan_vector<double>(size);
    const auto msols = init_nan_vector<double>(size);
    const auto xis = init_nan_vector<uint64_t>(size);
    const auto sdIds = sdIds_for_uninitialised_superdrops(size);

    const auto nandata = InitSupersData{
        initdata.solutes, sdgbxindexes, coord3s, coord1s, coord2s, radii, msols, xis, sdIds};

    return initdata + nandata;
  }

  /* sets sdIds for un-initialised superdrops' using an sdId's generator */
  std::vector<Superdrop::IDType> sdIds_for_uninitialised_superdrops(const size_t size) const {
    auto sdIdgen = Superdrop::IDType::Gen();

    auto sdIds = std::vector<Superdrop::IDType>(
        size, sdIdgen.set(std::numeric_limits<unsigned int>::signaling_NaN()));

    return sdIds;
  }

 public:
  /* constructor ensures the number of super-droplets to intialise is >= maxiumum number of
   * superdrops*/
  explicit InitSupersFromBinary(const OptionalConfigParams::InitSupersFromBinaryParams& config,
                                const GbxMaps& gbxmaps)
      : maxnsupers(config.maxnsupers),
        initnsupers(config.initnsupers),
        initsupers_filename(config.initsupers_filename),
        nspacedims(config.nspacedims),
        gbxmaps(gbxmaps) {
    if (maxnsupers < initnsupers) {
      const std::string err("cannot initialise more than the total number of super-droplets, ie. " +
                            std::to_string(maxnsupers) + " < " + std::to_string(initnsupers));
      throw std::invalid_argument(err);
    }
  }

  auto get_maxnsupers() const { return maxnsupers; }

  auto get_nspacedims() const { return nspacedims; }

  /* return InitSupersData created by reading data from a binary file to initialise "initnsupers"
  superdrops and then fills the rest of "maxnsupers" with un-initialised (and out of bounds)
  super-droplets. Also checks that the data created has the expected sizes. */
  InitSupersData fetch_data() const {
    auto initdata = fetch_superdrops_from_file();
    initdata = add_uninitialised_superdrops_data(initdata);
    trim_nonlocal_superdrops(initdata);
    check_initdata_sizes(initdata, maxnsupers, nspacedims);
    return initdata;
  }
};

#endif  // LIBS_INITIALISE_INIT_SUPERS_FROM_BINARY_HPP_
