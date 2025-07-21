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
 * Kokkos::parallel_for([...]) (on host) is equivalent to:
 * for (size_t ii(0); ii < ngbxs; ++ii){[...]} when in serial
 *
 */
void is_gbxinit_complete(const size_t ngbxs_from_maps, dualview_gbx gbxs,
                         const viewd_constsupers totsupers) {
  const auto ngbxs = size_t{gbxs.extent(0)};
  if (ngbxs != ngbxs_from_maps) {
    const std::string err(
        "number of gridboxes created not "
        "consistent with gridbox maps ie. " +
        std::to_string(ngbxs) + " != " + std::to_string(ngbxs_from_maps));
    throw std::invalid_argument(err);
  }

  const auto d_gbxs = gbxs.view_device();
  Kokkos::parallel_for(
      "is_gbxinit_complete", TeamPolicy(ngbxs, KCS::team_size),
      KOKKOS_LAMBDA(const TeamMember &team_member) {
        const auto ii = team_member.league_rank();
        assert(d_gbxs(ii).supersingbx.iscorrect(team_member, totsupers) &&
               "incorrect references to superdrops in gridbox");
      });
}

/**
 * @brief Print some information about initial Gridboxes.
 *
 * This function prints information about each Gridbox, including its index,
 * volume, and number of super-droplets.
 *
 */
void print_gbxs(const viewh_constgbx h_gbxs) {
  const auto ngbxs = size_t{h_gbxs.extent(0)};
  for (size_t ii(0); ii < ngbxs; ++ii) {
    const auto nsupers = h_gbxs(ii).supersingbx.nsupers();
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
 */
State GenGridbox::state_at(const unsigned int ii, const double volume) const {
  /* Type cast from std::pair to Kokkos::pair */
  const auto wvel = Kokkos::pair<double, double>{wvels.at(ii)};
  const auto uvel = Kokkos::pair<double, double>{uvels.at(ii)};
  const auto vvel = Kokkos::pair<double, double>{vvels.at(ii)};

  /* Return ii'th State from the initial conditions */
  return State(volume, presss.at(ii), temps.at(ii), qvaps.at(ii), qconds.at(ii), wvel, uvel, vvel);
}
