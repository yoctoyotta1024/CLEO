/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: streamout_observer.cpp
 * Project: observers
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Struct satisfies observer type and streams out live data to an output device
 * (e.g. computer screen) about the state of gridboxes during every observation
 * at fixed 'interval' timesteps.
 */

#include "observers/streamout_observer.hpp"

/**
 * @brief Prints a statement about the state of grid boxes.
 *
 * This function prints out information about the state of gridboxes.
 * It extracts information from the 0th gridbox in the gridboxes' view and prints some information
 * e.g. temperature, pressure, specific humidity, and specific cloud water content.
 * Additionally, it prints the total number of superdroplets in the domain and the total number of
 * gridboxes.
 *
 */
void StreamOutObserver::streamout_statement(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                                            const subviewd_constsupers d_supers) const {
  /* copy first gridbox into mirror view in case Gridboxes view is in device memory */
  auto d_gbx = Kokkos::subview(d_gbxs, kkpair_size_t({0, 1}));
  auto h_gbx = Kokkos::create_mirror_view_and_copy(HostSpace(), d_gbx);
  const auto &gbx = h_gbx(0);

  std::cout << "t=" << std::fixed << std::setprecision(2) << step2realtime(t_mdl)
            << "s, totnsupers=" << d_supers.extent(0) << ", ngbxs=" << d_gbxs.extent(0) << ", (Gbx"
            << gbx.get_gbxindex() << ": [T, p, qv, qc] = [" << gbx.state.temp * dlc::TEMP0 << "K, "
            << gbx.state.press * dlc::P0 << "Pa, " << std::scientific << std::setprecision(4)
            << gbx.state.qvap << ", " << gbx.state.qcond
            << "], nsupers = " << gbx.supersingbx.nsupers() << ")" << std::endl;
}
