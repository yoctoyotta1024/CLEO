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

void CreateGbxs::ensure_initialisation_complete(dualview_gbx gbxs,
                                                const size_t size) const
{
  const size_t ngbxs(gbxs.extent(0));

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
    if (!(gbxs.view_host()(ii).get_gbxindex() ==
          gbxs.view_device()(ii).get_gbxindex()))
    {
      const std::string err("gridboxes on device don't match host");
      throw std::invalid_argument(err);
    }
  
    if (!(gbxs.view_host()(ii).supersingbx.iscorrect()))
    {
      const std::string err("incorrect references to "
                            "superdrops in gridbox");
      throw std::invalid_argument(err);
    }
  }

}

void CreateGbxs::print_gbxs(dualview_gbx gbxs) const
/* print gridboxes information */
{
  for (size_t ii(0); ii < gbxs.extent(0); ++ii)
  {
    auto h_gbx(gbxs.view_host()(ii));
    const size_t nsupers(h_gbx.supersingbx.nsupers());
    std::cout << "gbx: " << h_gbx.get_gbxindex() 
    << ", vol = " << h_gbx.state.get_volume()
    << ", nsupers = " << nsupers << " \n";
    std::cout << "gbx: " << gbxs.view_device()(ii).get_gbxindex() << "\n";
  }
}