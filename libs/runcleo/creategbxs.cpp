/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: creategbxs.cpp
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
 * non-templated functionality for creating initialised view of Gridboxes from some
 * initial conditions
 */

#include "runcleo/creategbxs.hpp"

/**
 * @brief Check if gridbox initialisation is complete.
 *
 * This function checks if the number of created gridboxes is consistent with
 * the number of gridboxes from gridbox maps and if each gridbox has correct
 * references to superdroplets.
 *
 * @param ngbxs_from_maps The number of gridboxes from gridbox maps.
 * @param gbxs The dualview containing gridboxes.
 */
void is_gbxinit_complete(const size_t ngbxs_from_maps, dualview_gbx gbxs) {
  gbxs.sync_host();  // copy device to host (if sync flag was modified prior)
  const size_t ngbxs(gbxs.extent(0));
  const auto h_gbxs(gbxs.view_host());

  if (ngbxs != ngbxs_from_maps) {
    const std::string err(
        "number of gridboxes created not "
        "consistent with gridbox maps ie. " +
        std::to_string(ngbxs) + " != " + std::to_string(ngbxs_from_maps));
    throw std::invalid_argument(err);
  }

  for (size_t ii(0); ii < ngbxs; ++ii) {
    if (!(h_gbxs(ii).supersingbx.iscorrect())) {
      const std::string err(
          "incorrect references to "
          "superdrops in gridbox");
      throw std::invalid_argument(err);
    }
  }
}

/**
 * @brief Print some information about initial Gridboxes.
 *
 * This function prints information about each Gridbox, including its index,
 * volume, and number of super-droplets.
 *
 * @param h_gbxs The host view of Gridboxes.
 */
void print_gbxs(const viewh_constgbx h_gbxs) {
  const size_t ngbxs(h_gbxs.extent(0));
  for (size_t ii(0); ii < ngbxs; ++ii) {
    const size_t nsupers(h_gbxs(ii).supersingbx.nsupers());
    std::cout << "gbx: " << h_gbxs(ii).get_gbxindex()
              << ", (vol = " << h_gbxs(ii).state.get_volume() << ", nsupers = " << nsupers << ")\n";
  }
}

/**
 * @brief Get the state of a specified Gridbox from the initial conditions.
 *
 * This function returns the state of the Gridbox at the ii'th index in the initial conditions
 * given by the GenGridbox struct.
 *
 * @param ii The index of the Gridbox in the initial conditions data.
 * @param volume The volume of the Gridbox.
 * @return The State of the Gridbox.
 */
State GenGridbox::state_at(const unsigned int ii, const double volume) const {
  /* Type cast from std::pair to Kokkos::pair */
  Kokkos::pair<double, double> wvel(wvels.at(ii));
  Kokkos::pair<double, double> uvel(uvels.at(ii));
  Kokkos::pair<double, double> vvel(vvels.at(ii));

  /* Return ii'th State from the initial conditions */
  return State(volume, presss.at(ii), temps.at(ii), qvaps.at(ii), qconds.at(ii), wvel, uvel, vvel);
}
