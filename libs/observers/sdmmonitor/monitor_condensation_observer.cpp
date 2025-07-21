/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: monitor_condensation_observer.cpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to monitor condensation SDM microphysical process
 */

#include "./monitor_condensation_observer.hpp"

/**
 * @brief Parallel loop to fill d_data with zero value.
 */
void MonitorCondensation::reset_monitor() const {
  Kokkos::parallel_for(
      "reset_monitor", Kokkos::RangePolicy(0, d_data.extent(0)),
      KOKKOS_CLASS_LAMBDA(const size_t& jj) { d_data(jj) = 0.0; });
}

/**
 * @brief Monitor mass of liquid change due to condensation / evaporation
 *
 * Add totmass_condensed to current value for mass condensed since d_data was last reset.
 *
 * _Note:_ possible conversion of mass condensed at one timestep from double precision
 * (8 bytes double) to single precision (4 bytes float) in output depending on datatype alias.
 *
 * @param team_member Kokkkos team member in TeamPolicy parallel loop over gridboxes
 * @param totmass_condensed Mass condensed in one gridbox during one microphysical timestep
 */
KOKKOS_FUNCTION
void MonitorCondensation::monitor_condensation(const TeamMember& team_member,
                                               const double totmass_condensed) const {
  Kokkos::single(Kokkos::PerTeam(team_member), [=, this]() {
    const auto ii = team_member.league_rank();
    const auto mass_cond = static_cast<datatype>(totmass_condensed);
    d_data(ii) += mass_cond;
  });
}
