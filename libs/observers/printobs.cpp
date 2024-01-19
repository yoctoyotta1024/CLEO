/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: printobs.cpp
 * Project: observers
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer Concept and related structures for various ways
 * of observing (outputing data from) CLEO.
 * An example of an observer is printing some data
 * from a gridbox's thermostate to the terminal
 */


#include "./printobs.hpp"

void PrintObserver::
    print_statement(const unsigned int t_mdl,
                    const viewh_constgbx h_gbxs,
                    const viewd_constsupers totsupers) const
{
  const auto gbx = h_gbxs(0);
  std::cout << "t="
            << std::fixed << std::setprecision(2)
            << step2realtime(t_mdl)
            << "s, totnsupers=" << totsupers.extent(0)
            << ", ngbxs=" << h_gbxs.extent(0)
            << ", (Gbx" << gbx.get_gbxindex()
            << ": [T, p, qv, qc] = [" << gbx.state.temp * dlc::TEMP0
            << "K, " << gbx.state.press * dlc::P0 << "Pa, "
            << std::scientific << std::setprecision(4)
            << gbx.state.qvap
            << ", " << gbx.state.qcond
            << "], nsupers = " << gbx.supersingbx.nsupers() << ")\n";
}
