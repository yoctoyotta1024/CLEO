/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: gensuperdrop.cpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality to generate a superdroplet (on device) from some initial conditions.
 */

#include "runcleo/gensuperdrop.hpp"

/**
 * @brief Generate a super-droplet using initial data for the kk'th superdrop.
 *
 * This function returns a superdrop generated from the specified position
 * 'kk' in the initial conditions data.
 *
 * @param kk The index of the superdrop to generate.
 * @return The generated super-droplet.
 */
Superdrop GenSuperdrop::operator()(const unsigned int kk) const {
  const auto sdgbxindex = initdata.sdgbxindexes.at(kk);
  const auto coords312 = coords_at(kk);
  const SuperdropAttrs attrs(attrs_at(kk));
  const auto sdId = initdata.sdIds.at(kk);

  return Superdrop(sdgbxindex, coords312[0], coords312[1], coords312[2], attrs, sdId);
}

/**
 * @brief Returns initial spatial coordinates of the kk'th super-droplet.
 *
 * A coordinate is only copied from the corresponding coords vector if that coordinate is
 * consistent with the number of spatial dimensions of the model. Otherwise, the coordinate is set
 * to 0.0. For example, if the model is 1-D, only coord3 is obtained from the initial data vector;
 * coord1 and coord2 are set to 0.0.
 *
 * @param kk The index of the super-droplet in the initial
 * data (0 <= kk < total number of super-droplets).
 * @return An array containing the kk'th superdroplet's spatial
 * coordinates (coord3, coord1, coord2).
 */
std::array<double, 3> GenSuperdrop::coords_at(const unsigned int kk) const {
  std::array<double, 3> coords312{0.0, 0.0, 0.0};

  switch (nspacedims) {
    case 3:  // 3-D model
      coords312[2] = initdata.coord2s.at(kk);
      [[fallthrough]];
    case 2:  // 3-D or 2-D model
      coords312[1] = initdata.coord1s.at(kk);
      [[fallthrough]];
    case 1:  // 3-D, 2-D or 1-D model
      coords312[0] = initdata.coord3s.at(kk);
  }

  return coords312;
}

/* helper function to return a superdroplet's attributes
at position kk in the initial conditions data. All
superdroplets created with same solute properties */
/**
 * @brief Function returns a superdroplet's attributes
 * from position 'kk' in the initial conditions data. All
 * super-droplets have the same solute properties.
 *
 * @param kk The index of the super-droplet in the initial
 * data (0 <= kk < total number of super-droplets).
 * @return The attributes of the superdrop from index 'kk'.
 */
SuperdropAttrs GenSuperdrop::attrs_at(const unsigned int kk) const {
  const auto radius = initdata.radii.at(kk);
  const auto msol = initdata.msols.at(kk);
  const auto xi = initdata.xis.at(kk);
  const SoluteProperties solute(initdata.solutes.at(0));

  return SuperdropAttrs(solute, xi, radius, msol, true);
}
