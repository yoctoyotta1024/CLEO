/*
 * ----- CLEO -----
 * File: creategbxs.hpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 9th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * file for structure to create a dualview of
 * gridboxes from using some initial conditions
 */

#ifndef CREATEGBXS_HPP
#define CREATEGBXS_HPP

#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <stdexcept>
#include <utility>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_DualView.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gbxindex.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/supersingbx.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/state.hpp"

template <GridboxMaps GbxMaps, typename GbxInitConds>
dualview_gbx create_gbxs(const GbxMaps &gbxmaps,
                         const GbxInitConds &gbxic,
                         const viewd_supers totsupers);

class GenGridbox
{
private:
  std::shared_ptr<Gbxindex::Gen> GbxindexGen; // pointer to gridbox index generator
  std::vector<double> presss;
  std::vector<double> temps;
  std::vector<double> qvaps;
  std::vector<double> qconds;
  std::vector<std::pair<double, double>> wvels;
  std::vector<std::pair<double, double>> uvels;
  std::vector<std::pair<double, double>> vvels;

  KOKKOS_FUNCTION
  State state_at(const unsigned int ii,
                 const double volume) const;

public:
  template <typename GbxInitConds>
  GenGridbox(const GbxInitConds &gbxic)
      : GbxindexGen(std::make_shared<Gbxindex::Gen>()),
        presss(gbxic.press()),
        temps(gbxic.temp()),
        qvaps(gbxic.qvap()),
        qconds(gbxic.qcond()),
        wvels(gbxic.wvel()),
        uvels(gbxic.uvel()),
        vvels(gbxic.vvel()) {}

  template <GridboxMaps GbxMaps>
  Gridbox operator()(const unsigned int ii,
                     const GbxMaps &gbxmaps,
                     const viewd_supers totsupers) const
  {
    const auto gbxindex(GbxindexGen->next(ii));
    const double volume(gbxmaps.get_gbxvolume(gbxindex.value));
    const State state(state_at(ii, volume));

    return Gridbox(gbxindex, state, totsupers);
  }

  template <GridboxMaps GbxMaps>
  Gridbox operator()(const unsigned int ii,
                     const GbxMaps &gbxmaps,
                     const viewd_supers totsupers,
                     const SupersInGbx::kkpair refs) const
  {
    const auto gbxindex(GbxindexGen->next(ii));
    const double volume(gbxmaps.get_gbxvolume(gbxindex.value));
    const State state(state_at(ii, volume));
    
    return Gridbox(gbxindex, state, totsupers, refs);
  }
};

template <GridboxMaps GbxMaps>
inline void initialise_gbxs_on_host(const GbxMaps &gbxmaps,
                                    const GenGridbox &gen,
                                    const viewd_supers totsupers,
                                    const viewh_gbx h_gbxs);
/* initialise the host view of gridboxes
using some data from a GbxInitConds instance
e.g. for each gridbox's volume */

template <GridboxMaps GbxMaps, typename GbxInitConds>
inline dualview_gbx initialise_gbxs(const GbxMaps &gbxmaps,
                                    const GbxInitConds &gbxic,
                                    const viewd_supers totsupers);
/* initialise a dualview of gridboxes (on host and device
memory) using data from a GbxInitConds instance to initialise
the host view and then syncing the view to the device */

void is_gbxinit_complete(const size_t ngbxs_from_maps,
                         dualview_gbx gbxs);

void print_gbxs(const viewh_constgbx gbxs);
/* print gridboxes information */

template <GridboxMaps GbxMaps, typename GbxInitConds>
dualview_gbx create_gbxs(const GbxMaps &gbxmaps,
                         const GbxInitConds &gbxic,
                         const viewd_supers totsupers)
{

  std::cout << "\n--- create gridboxes ---\n"
            << "initialising\n";
  const dualview_gbx gbxs(initialise_gbxs(gbxmaps, gbxic, totsupers));

  std::cout << "checking initialisation\n";
  is_gbxinit_complete(gbxmaps.maps_size() - 1, gbxs);
  print_gbxs(gbxs.view_host());

  std::cout << "--- create gridboxes: success ---\n";

  return gbxs;
}

template <GridboxMaps GbxMaps, typename GbxInitConds>
inline dualview_gbx
initialise_gbxs(const GbxMaps &gbxmaps,
                const GbxInitConds &gbxic,
                const viewd_supers totsupers)
/* initialise a view of superdrops (on device memory)
using data from an InitData instance for their initial
gbxindex, spatial coordinates and attributes */
{
  // create dualview for gridboxes on device and host memory
  dualview_gbx gbxs("gbxs", gbxic.get_ngbxs());

  // initialise gridboxes on host
  const GenGridbox gen(gbxic);
  gbxs.sync_host();
  initialise_gbxs_on_host(gbxmaps, gen, totsupers, gbxs.view_host());
  gbxs.modify_host();

  // update device gridbox view to match host's gridbox view
  gbxs.sync_device();

  return gbxs;
}

template <GridboxMaps GbxMaps>
inline void initialise_gbxs_on_host(const GbxMaps &gbxmaps,
                                    const GenGridbox &gen,
                                    const viewd_supers totsupers,
                                    const viewh_gbx h_gbxs)
/* initialise the host view of gridboxes using
some data from a GbxInitConds instance e.g. for
each gridbox's volume. 
Kokkos::parallel_for([...]) is equivalent to:
for (size_t ii(0); ii < ngbxs; ++ii) {[...]}
when in serial */
{
  const size_t ngbxs(h_gbxs.extent(0));
  
  /* equivalent serial version of parallel_for loop below
  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    h_gbxs(ii) = gen(ii, gbxmaps, totsupers);
  }
  */

  Kokkos::parallel_for(
      "initialise_gbxs_on_host",
      HostTeamPolicy(ngbxs, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const HostTeamMember &team_member) {
        const int ii = team_member.league_rank();
        
        const SupersInGbx::kkpair refs = ?
        const Gridbox gbx(gen(ii, gbxmaps, totsupers, refs));

        Kokkos::single(
            Kokkos::PerTeam(team_member),
            [=]()
            {
              h_gbxs(ii) = gbx;
            });
      });
}

#endif // CREATEGBXS_HPP