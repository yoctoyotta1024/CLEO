/*
 * ----- CLEO -----
 * File: initsupers_frombinary.hpp
 * Project: initialise
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 2nd November 2023
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

#ifndef INITSUPERS_FROMBINARY_HPP
#define INITSUPERS_FROMBINARY_HPP

#include <vector>
#include <string_view>

#include "./config.hpp"
#include "./readbinary.hpp"
#include "superdrops/superdrop_attrs.hpp"

struct InitSupersFromBinary
/* struct containing functions which return data
for the initial conditions needed to create
superdroplets e.g. via the CreateSupers struct */
{
private:
  size_t totnsupers; // total number of superdroplets (in kokkos view on device initially)
  unsigned int nspacedims; // number of spatial dimensions to model (0-D, 1-D, 2-D of 3-D)
  std::string_view initsupers_filename; // name of binary file for some of superdrops' initial conditons

  void init_solutes_data(InitSupersData &initdata) const;
  /* sets initial data for solutes as
  a single SoluteProprties instance */

  void initdata_from_binary(InitSupersData &initdata) const;
  /* sets initial data in initdata using data read
  from a binary file called initsupers_filename */

  void read_initdata_binary(InitSupersData &initdata,
                            std::ifstream &file,
                            const std::vector<VarMetadata> &meta) const;
  /* copy data for vectors from binary file to initdata struct */

  void check_initdata_sizes(const InitSupersData &initdata) const;
  /* check all the vectors in the initdata struct all
  have sizes consistent with one another. Include
  coords data in check if nspacedims != 0 */    

public:
  InitSupersFromBinary(const Config &config)
      : totnsupers(config.totnsupers),
        nspacedims(config.nspacedims),
        initsupers_filename(config.initsupers_filename) {}

  auto get_totnsupers() const { return totnsupers; }

  auto get_nspacedims() const { return nspacedims; }

  size_t fetch_data_size() const; 
  /* data size returned is number of variables as
  declared by the metadata for the first variable
  in the initsupers file */

  InitSupersData fetch_data() const;
  /* return InitSupersData created by reading a binary
  file and creating a SoluteProperties struct.
  Then check that the input data has the correct sizes. */
};

#endif // INITSUPERS_FROMBINARY_HPP