/*
 * ----- CLEO -----
 * File: creategbxs.cpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 */

#include "./creategbxs.hpp"

void CreateGbxs::ensure_initialisation_complete(dualview_gbx supers,
                                                const size_t size) const
{
}

void CreateGbxs::print_gbxs(dualview_gbx gbxs) const
/* print gridboxes information */
{
  for (size_t ii(0); ii < gbxs.extent(0); ++ii)
  {
    std::cout << "gbx: " << gbxs.view_host()(ii).get_gbxindex() << "\n";
    std::cout << "gbx: " << gbxs.view_device()(ii).get_gbxindex() << "\n";
  }
}