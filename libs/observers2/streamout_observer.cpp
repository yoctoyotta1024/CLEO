/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: streamout_observer.cpp
 * Project: observers2
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 27th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Struct satisfies observer type and streams out live data to an output device
 * (e.g. computer screen) about the state of gridboxes during every observation
 * at fixed 'interval' timesteps.
 */

#include "observers2/streamout_observer.hpp"

void StreamOutObserver::print_statement(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                                        const viewd_constsupers totsupers) const {
  /* copy first gridbox into mirror view in case Gridboxes view is in device memory */
  auto h_gbx = Kokkos::create_mirror_view(Kokkos::subview(d_gbxs, kkpair_size_t({0, 1})));
  Kokkos::deep_copy(h_gbx, d_gbxs(0));
  const auto &gbx(h_gbx(0));

  std::cout << "t=" << std::fixed << std::setprecision(2) << step2realtime(t_mdl)
            << "s, totnsupers=" << totsupers.extent(0) << ", ngbxs=" << d_gbxs.extent(0) << ", (Gbx"
            << gbx.get_gbxindex() << ": [T, p, qv, qc] = [" << gbx.state.temp * dlc::TEMP0 << "K, "
            << gbx.state.press * dlc::P0 << "Pa, " << std::scientific << std::setprecision(4)
            << gbx.state.qvap << ", " << gbx.state.qcond
            << "], nsupers = " << gbx.supersingbx.nsupers() << ")\n";
}
