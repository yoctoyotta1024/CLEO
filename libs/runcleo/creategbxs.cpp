/*
 * ----- CLEO -----
 * File: creategbxs.cpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 29th October 2023
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

void is_gbxinit_complete(dualview_gbx gbxs,
                             const size_t size)
{
  gbxs.sync_host(); // copy device to host (if prior flag was set)
  const size_t ngbxs(gbxs.extent(0));
  const auto h_gbxs(gbxs.view_host());
 
  if (!(ngbxs == size))
  {
    const std::string err("number of gridboxes created not "
                          "consistent with initialisation data ie. " +
                          std::to_string(ngbxs) + " != " +
                          std::to_string(size));
    throw std::invalid_argument(err);
  }

  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    if (!(h_gbxs(ii).supersingbx.iscorrect()))
    {
      const std::string err("incorrect references to "
                            "superdrops in gridbox");
      throw std::invalid_argument(err);
    }
  }
}

void print_gbxs(const viewh_constgbx h_gbxs)
/* print gridboxes information */
{
  const size_t ngbxs(h_gbxs.extent(0));
  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    const size_t nsupers(h_gbxs(ii).supersingbx.nsupers());
    std::cout << "gbx: " << h_gbxs(ii).get_gbxindex()
              << ", (vol = " << h_gbxs(ii).state.get_volume()
              << ", nsupers = " << nsupers << ")\n";
  }
}