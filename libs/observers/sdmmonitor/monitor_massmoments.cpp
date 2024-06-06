/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_massmoments.cpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 6th June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to monitor mass moments from SDM microphysical process
 */

#include "./monitor_massmoments.hpp"

/**
 * @brief Parallel loop to fill d_data with zero value
 */
void MonitorMassMomentViews::reset_views() const {
  Kokkos::parallel_for(
      "reset_views", Kokkos::RangePolicy(0, d_mom0.extent(0)),
      KOKKOS_CLASS_LAMBDA(const size_t& jj) {
        d_mom0(jj) = 0;
        d_mom1(jj) = 0.0;
        d_mom2(jj) = 0.0;
      });
}

/**
 * @brief Write the 0th, 1st and 2nd moments of the droplet mass distribution to data views.
 *
 * Calculates the current mass moments and then overwrites the current values for the
 * mass moments stored since the data views were last reset.
 *
 * _Note:_ possible conversion of mass moments at one timestep from double precision
 * (8 bytes double) to single precision (4 bytes float) in output.
 *
 * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
 * @param supers (sub)View of all the superdrops in one gridbox
 */
KOKKOS_FUNCTION
void MonitorMassMomentViews::calculate_massmoments(const TeamMember& team_member,
                                                   const viewd_constsupers supers) const {
  const auto ii = team_member.league_rank();  // position of gridbox

  d_mom0(ii) =
      static_cast<uint64_t>(ii);  // TODO(CB): WIP call observer.hpp calculate_massmoments instead
  d_mom1(ii) = static_cast<float>(10.0 * ii);
  d_mom2(ii) = static_cast<float>(100.0 * ii);
}
