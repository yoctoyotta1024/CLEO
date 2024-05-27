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
 * Last Modified: Monday 27th May 2024
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
 * @brief Parallel loop to fill d_data with zero value and set monitor_count to zero.
 */
void MonitorMassMoments::reset_monitor() const {
  Kokkos::parallel_for(
      "reset_monitor", Kokkos::RangePolicy(0, d_data.extent(0)),
      KOKKOS_CLASS_LAMBDA(const size_t& jj) { d_data(jj) = 0.0; });
  monitor_count = 0;
}

/**
 * @brief Monitor 0th, 1st and 2nd moments of the droplet mass distribution
 *
 * Averages current of mass moments with value stored since d_data was last reset.
 *
 * _Note:_ possible conversion of mass moments at one timestep from double precision
 * (8 bytes double) to single precision (4 bytes float) in output depending on datatype alias.
 *
 * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
 * @param supers (sub)View of all the superdrops in one gridbox during one microphysical timestep
 */
KOKKOS_FUNCTION
void MonitorMassMoments::monitor_microphysics(const TeamMember& team_member,
                                              const viewd_constsupers supers) const {
  const auto ii = team_member.league_rank();  // position of gridbox
  const auto massmom = ii / monitor_count;    // TODO(CB): calc mass moments properly
  d_data(ii) += massmom;

  ++monitor_count;
}
