/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: init_all_supers_from_binary.hpp
 * Project: initialise
 * Created Date: Tuesday 17th October 2023
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

#ifndef LIBS_INITIALISE_INIT_ALL_SUPERS_FROM_BINARY_HPP_
#define LIBS_INITIALISE_INIT_ALL_SUPERS_FROM_BINARY_HPP_

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "configuration/optional_config_params.hpp"
#include "initialise/initialconditions.hpp"
#include "initialise/readbinary.hpp"
#include "superdrops/superdrop.hpp"

/* check all the vectors in the initdata struct all have sizes consistent with one another
and with maxnsupers. Include coords data in check if nspacedims > 0 */
void check_initdata_sizes(const InitSupersData &in, const size_t maxnsupers,
                          const size_t nspacedims);

/* struct containing functions which return data
for the initial conditions needed to create
superdroplets e.g. via the CreateSupers struct */
struct InitAllSupersFromBinary {
 private:
  size_t maxnsupers; /**< total number of super-droplets (in kokkos view on device) */
  std::filesystem::path initsupers_filename; /**< filename for super-droplets' initial conditons */
  unsigned int nspacedims; /**< number of spatial dimensions to model (0-D, 1-D, 2-D of 3-D) */

  /* sets initial data for solutes as
  a single SoluteProprties instance */
  void initdata_for_solutes(InitSupersData &initdata) const;

  /* sets initial data for sdIds using its generator */
  void initdata_for_sdIds(InitSupersData &initdata) const;

  /* sets initial data in initdata using data read
  from a binary file called initsupers_filename */
  void initdata_from_binary(InitSupersData &initdata) const;

  /* copy data for vectors from binary file to initdata struct */
  void read_initdata_binary(InitSupersData &initdata, std::ifstream &file,
                            const std::vector<VarMetadata> &meta) const;

  /* data size returned is number of variables as
  declared by the metadata for the first variable
  in the initsupers file */
  size_t fetch_data_size() const;

 public:
  /**
   * @brief Constructor for InitAllSupersFromBinary.
   *
   * Calls constructor for InitAllSupersFromBinary wiht additional assert to sanity check
   * condfiguration matches initialisation expected by this struct.
   *
   * @param config Configuration for member variables.
   *
   */
  explicit InitAllSupersFromBinary(const OptionalConfigParams::InitSupersFromBinaryParams &config)
      : InitAllSupersFromBinary(config.maxnsupers, config.initsupers_filename, config.nspacedims) {
    assert((config.maxnsupers == config.initnsupers) &&
           "configuration parameter not consistent with "
           "initialising all super-droplets from binary");
  }

  /**
   * @brief Constructor for InitAllSupersFromBinary.
   *
   * This function checks if there is enough data in the initialisation files for
   * super-droplets in order to initialise "maxnsupers" superdroplets. If the initialisation
   * data is the wrong size, it throws an exception with the appropriate error message.
   *
   * @param maxnsupers The expected initial total number of super-droplets
   * @param initsupers_filename filename for super-droplets' initial conditons
   * @param nspacedims Number of spatial dimensions to model (0-D, 1-D, 2-D of 3-D)
   *
   * @throws std::invalid_argument If the number of super-droplets is wrong.
   */
  InitAllSupersFromBinary(const size_t maxnsupers, const std::filesystem::path initsupers_filename,
                          const unsigned int nspacedims)
      : maxnsupers(maxnsupers), initsupers_filename(initsupers_filename), nspacedims(nspacedims) {
    const auto size = fetch_data_size();

    if (maxnsupers < size) {
      const std::string err(
          "Fewer superdroplets will be created than is given by initialisation data, ie. " +
          std::to_string(maxnsupers) + " < " + std::to_string(size));
      throw std::invalid_argument(err);
    } else if (maxnsupers > size) {
      const std::string err(
          "Not enough initialisation data for number of superdroplets that will be created, ie." +
          std::to_string(maxnsupers) + " > " + std::to_string(size));
      throw std::invalid_argument(err);
    }

    assert((maxnsupers == size) && "size of data not the same as the number of superdroplets");
  }

  auto get_maxnsupers() const { return maxnsupers; }

  auto get_nspacedims() const { return nspacedims; }

  /* returns InitSupersData created by reading a binary file and creating a
  SoluteProperties struct. Also checks that the data created has the expected sizes. */
  InitSupersData fetch_data() const {
    auto initdata = InitSupersData{};

    initdata_for_solutes(initdata);
    initdata_for_sdIds(initdata);
    initdata_from_binary(initdata);
    check_initdata_sizes(initdata, maxnsupers, nspacedims);

    return initdata;
  }
};

#endif  // LIBS_INITIALISE_INIT_ALL_SUPERS_FROM_BINARY_HPP_
