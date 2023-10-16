/*
 * ----- CLEO -----
 * File: runcleo.cpp
 * Project: runcleo
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 14th October 2023
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

dualview_gbx create_gridboxes()
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

viewd_supers create_superdrops()
/* create view of superdrops */
{
  const int nSDsvec(25);
  viewd_supers supers("SDs", nSDsvec);
  
  const auto solute(std::make_shared<const SoluteProperties>());
  auto sdIdGen = Superdrop::IDType::Gen{};
  for (int jj(0); jj < nSDsvec; ++jj)
  {
    const unsigned int sd_gbxindex = jj;
    const auto sd_id = sdIdGen.next();
    const SuperdropAttrs attrs{0.0, 0.0, 0, solute};
    const std::vector<double> crds{0.0, 0.0, 0.0};

    supers(jj) = Superdrop(sd_gbxindex, crds.at(0), crds.at(1),
                           crds.at(2), attrs, sd_id);
  }
  
  for (int jj(0); jj < nSDsvec; ++jj)
  { 
    std::cout << "sdid: " << supers(jj).id.value << "\n"; 
  }

  return supers;
}

unsigned int next_stepsize(const unsigned int t_mdl,
                           const SDMMethods &sdm)
/* returns size of next step of model ('onestep')
given current time t_mdl, so that next time
(t_next = t_mdl + onestep) is time of obs or coupl */
{
  const unsigned int couplstep(sdm.get_couplstep());
  const unsigned int obsstep(sdm.obs.get_obsstep());

  const auto next_step = [t_mdl](const unsigned int interval)
  {
    return ((t_mdl / interval) + 1) * interval;
  };

  /* t_next is smaller out of time of next coupl and obs */
  const unsigned int next_coupl(next_step(couplstep));
  const unsigned int next_obs(next_step(obsstep));

  return std::min(next_coupl, next_obs) - t_mdl;
}

