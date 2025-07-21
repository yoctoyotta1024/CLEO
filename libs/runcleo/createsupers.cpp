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
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * non-templated functionality required by RunCLEO to create a view of superdroplets (on device)
 * using some initial conditions
 */

#include "runcleo/createsupers.hpp"

/**
 * @brief Check if superdroplets initialisation is complete.
 *
 * This function checks if the initialisation of supers view is complete by checking if the
 * superdroplets are sorted by ascending gridbox indexes. If the initialisation is incomplete
 * (the superdroplets are not sorted), it throws an exception with an appropriate error message.
 *
 * @param allsupers Struct to handle all the super-droplets in device memory.
 *
 * @throws std::invalid_argument If the initialisation is incomplete i.e. the super-droplets
 * are not ordered correctly.
 */
void is_sdsinit_complete(const SupersInDomain allsupers) {
  if (allsupers.is_sorted() == 0) {
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
 * @param totsupers The view of super-droplets in device memory.
 */
void print_supers(const viewd_constsupers totsupers) {
  auto h_totsupers =
      Kokkos::create_mirror_view(totsupers);  // mirror of supers in case view is on device memory
  Kokkos::deep_copy(h_totsupers, totsupers);

  for (size_t kk(0); kk < h_totsupers.extent(0); ++kk) {
    std::cout << "SD: " << h_totsupers(kk).sdId << " [gbx, (coords), (attrs)]: [ "
              << h_totsupers(kk).get_sdgbxindex() << ", (" << h_totsupers(kk).get_coord3() << ", "
              << h_totsupers(kk).get_coord1() << ", " << h_totsupers(kk).get_coord2() << "), ("
              << h_totsupers(kk).is_solute() << ", " << h_totsupers(kk).get_radius() << ", "
              << h_totsupers(kk).get_msol() << ", " << h_totsupers(kk).get_xi() << ") ] \n";
  }
}
