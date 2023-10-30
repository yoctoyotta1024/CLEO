/*
 * ----- CLEO -----
 * File: initsupers_frombinary.hpp
 * Project: initialise
 * Created Date: Tuesday 17th October 2023
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

#ifndef INITSUPERS_FROMBINARY_HPP
#define INITSUPERS_FROMBINARY_HPP

#include <vector>
#include <string_view>

#include "./config.hpp"
#include "runcleo/initialconditions.hpp"
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

  void init_solutes(InitSupersData &initdata) const;

  void initdata_from_binary(InitSupersData &initdata) const; 

public:
  InitSupersFromBinary(const Config &config)
      : totnsupers(config.totnsupers),
        nspacedims(config.nspacedims),
        initsupers_filename(config.initsupers_filename) {}

  auto get_totnsupers() const { return totnsupers; }

  auto get_nspacedims() const { return nspacedims; }

  size_t fetch_data_size() const;

  InitSupersData fetch_data() const;
};

#endif // INITSUPERS_FROMBINARY_HPP