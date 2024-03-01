/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
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
 * File Description:
 * various structures used to make an
 * InitConds instance which satifies
 * the InitialConditions concept
 */

#ifndef LIBS_INITIALISE_INITCONDS_HPP_
#define LIBS_INITIALISE_INITCONDS_HPP_

#include <array>
#include <vector>

#include "superdrops/superdrop_attrs.hpp"

struct InitSupersData {
  std::array<SoluteProperties, 1> solutes;
  std::vector<unsigned int> sdgbxindexes;
  std::vector<double> coord3s;
  std::vector<double> coord1s;
  std::vector<double> coord2s;
  std::vector<double> radii;
  std::vector<double> msols;
  std::vector<uint64_t> xis;
};

/* struct for functions to generate
initial conditions for CLEO */
template <typename SuperdropInitConds, typename GbxInitConds>
struct InitConds {
  SuperdropInitConds initsupers;  // initial conditions for creating superdroplets
  GbxInitConds initgbxs;          // initial conditions for creating gridboxes

  InitConds(const SuperdropInitConds initsupers, const GbxInitConds initgbxs)
      : initsupers(initsupers), initgbxs(initgbxs) {}
};

#endif  // LIBS_INITIALISE_INITCONDS_HPP_
