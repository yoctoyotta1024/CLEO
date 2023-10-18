/*
 * ----- CLEO -----
 * File: runcleo.cpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
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
 * functionality related to timestepping CLEO coupled model
 * (CLEO SDM coupled one-way/two-ways to a Dynamics Solver)
 */


#include "./runcleo.hpp"

create_gridboxes()
/* create dualview of gridboxes (in general this
is two distinct views on host and device memory) */
{
  const size_t ngbxs(10);
  dualview_gbx gbxs("Gbxs", ngbxs);

  gbxs.sync_host();
  viewh_gbx h_gbxs = gbxs.view_host();

  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    const unsigned int gbxindex = ii;
    const double volume = 0.0;
    const Kokkos::pair<size_t, size_t> pos = {0, 3};
    
    h_gbxs(ii) = Gridbox(gbxindex, volume, pos);
  }
  gbxs.modify_host();
  gbxs.sync_device();

  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    std::cout << "gbx: " << gbxs.view_host()(ii).get_gbxindex() << "\n"; 
    std::cout << "gbx: " << gbxs.view_device()(ii).get_gbxindex() << "\n"; 
  }

  return gbxs;
}