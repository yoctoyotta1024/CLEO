/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: createsupers.cpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 8th February 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * non-templated functionality required by RunCLEO to create a view of superdroplets (on device)
 * using some initial conditions
 */

#include "runcleo/createsupers.hpp"

/**
 * @brief Check if superdroplets initialisation is complete and sorted.
 *
 * This function checks if the initialisation of supers view is complete and if the super-droplets
 * are sorted by ascending superdroplet Gridbox index. If the initialisation is incomplete
 * or the superdroplets are not sorted, it throws an exception with an appropriate error message.
 *
 * @param supers The view of super-droplets in device memory.
 * @param size The expected number of super-droplets.
 *
 * @throws std::invalid_argument If the initialisation is incomplete or super-droplets
 * are not ordered correctly.
 */
void is_sdsinit_complete(const viewd_constsupers supers, const size_t size) {
  if (supers.extent(0) < size) {
    const std::string err(
        "Fewer superdroplets were created than were"
        " given by initialisation data ie. " +
        std::to_string(supers.extent(0)) + " < " + std::to_string(size));
    throw std::invalid_argument(err);
  }

  if (is_sorted(supers) == 0) {
    const std::string err(
        "supers ordered incorrectly "
        "(ie. not sorted by asceding sdgbxindex");
    throw std::invalid_argument(err);
  }
}

/**
 * @brief Print statement about initialised super-droplets.
 *
 * This function prints information about each superdroplet, including its ID, Gridbox index,
 * spatial coordinates, and attributes.
 *
 * @param supers The view of super-droplets in device memory.
 */
void print_supers(const viewd_constsupers supers) {
  auto h_supers =
      Kokkos::create_mirror_view(supers);  // mirror of supers in case view is on device memory
  Kokkos::deep_copy(h_supers, supers);

  for (size_t kk(0); kk < h_supers.extent(0); ++kk) {
    std::cout << "SD: " << h_supers(kk).sdId.value << " [gbx, (coords), (attrs)]: [ "
              << h_supers(kk).get_sdgbxindex() << ", (" << h_supers(kk).get_coord3() << ", "
              << h_supers(kk).get_coord1() << ", " << h_supers(kk).get_coord2() << "), ("
              << h_supers(kk).is_solute() << ", " << h_supers(kk).get_radius() << ", "
              << h_supers(kk).get_msol() << ", " << h_supers(kk).get_xi() << ") ] \n";
  }
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
 * data (0 <= kk < total number of superdrops).
 * @return An array containing the kk'th superdroplet's spatial
 * coordinates (coord3, coord1, coord2).
 */
std::array<double, 3> GenSuperdrop::coords_at(const unsigned int kk) const {
  std::array<double, 3> coords312{0.0, 0.0, 0.0};

  switch (nspacedims) {
    case 3:  // 3-D model
      coords312[2] = initdata.coord2s.at(kk);
    case 2:  // 3-D or 2-D model
      coords312[1] = initdata.coord1s.at(kk);
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
 * data (0 <= kk < total number of superdrops).
 * @return The attributes of the superdrop from index 'kk'.
 */
SuperdropAttrs GenSuperdrop::attrs_at(const unsigned int kk) const {
  const auto radius = initdata.radii.at(kk);
  const auto msol = initdata.msols.at(kk);
  const auto xi = initdata.xis.at(kk);
  const SoluteProperties solute(initdata.solutes.at(0));

  return SuperdropAttrs(solute, xi, radius, msol);
}

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
  const auto sdId = sdIdGen->next(static_cast<size_t>(kk));

  return Superdrop(sdgbxindex, coords312[0], coords312[1], coords312[2], attrs, sdId);
}
