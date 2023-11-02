/*
 * ----- CLEO -----
 * File: initconds.hpp
 * Project: initialise
 * Created Date: Thursday 2nd November 2023
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
 */

#ifndef INITCONDS_HPP
#define INITCONDS_HPP

#include <array>
#include <vector>

#include "superdrops/superdrop_attrs.hpp"

struct InitSupersData
{
  std::array<SoluteProperties, 1> solutes;
  std::vector<unsigned int> sdgbxindexes;
  std::vector<double> coord3s;
  std::vector<double> coord1s;
  std::vector<double> coord2s;
  std::vector<double> radii;
  std::vector<double> msols;
  std::vector<unsigned long long> xis;
};

template <typename SuperdropInitConds, typename GbxInitConds>
struct InitConds
/* struct for functions to generate
initial conditions for CLEO */
{
  SuperdropInitConds initsupers; // initial conditions for creating superdroplets
  GbxInitConds initgbxs;         // initial conditions for creating gridboxes

  InitConds(const SuperdropInitConds initsupers,
            const GbxInitConds initgbxs)
      : initsupers(initsupers),
        initgbxs(initgbxs) {}
};

#endif // INITCONDS_HPP