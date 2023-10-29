/*
 * ----- CLEO -----
 * File: creategbxs.hpp
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

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_DualView.hpp>

#include "../kokkosaliases.hpp"
#include "./initconds.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/state.hpp"

template <typename GbxInitConds>
dualview_gbx create_gbxs(const GbxInitConds &gic,
                        const viewd_supers supers);

template <typename GbxInitConds>
inline void initialise_gbxs_on_host(const GbxInitConds &gic,
                                    const viewd_supers supers,
                                    const viewh_gbx h_gbxs);
/* initialise the host view of gridboxes
using some data from a GbxInitConds instance
e.g. for each gridbox's volume */

template <typename GbxInitConds>
inline dualview_gbx initialise_gbxs(const GbxInitConds &gic,
                                    const viewd_supers supers);
/* initialise a dualview of gridboxes (on host and device
memory) using data from a GbxInitConds instance to initialise
the host view and then syncing the view to the device */

void is_gbxinit_complete(dualview_gbx gbxs,
                                    const size_t size);

void print_gbxs(const viewh_constgbx gbxs);
/* print gridboxes information */

class GenGridbox
{
private:
  std::unique_ptr<Gridbox::Gbxindex::Gen> GbxindexGen; // pointer to gridbox index generator
  std::vector<double> volumes;
  std::vector<double> presss;
  std::vector<double> temps;
  std::vector<double> qvaps;
  std::vector<double> qconds;
  std::vector<Kokkos::pair<double, double>> wvels;
  std::vector<Kokkos::pair<double, double>> uvels;
  std::vector<Kokkos::pair<double, double>> vvels;

  State state_at(const unsigned int ii) const;

public:
  template <typename GbxInitConds>
  GenGridbox(const GbxInitConds &gic)
      : GbxindexGen(std::make_unique<Gridbox::Gbxindex::Gen>()),
        volumes(gic.volume()),
        presss(gic.press()),
        temps(gic.temp()),
        qvaps(gic.qvap()),
        qconds(gic.qcond()),
        wvels(gic.wvel()),
        uvels(gic.uvel()),
        vvels(gic.vvel()) {}

  Gridbox operator()(const unsigned int ii,
                      const viewd_supers supers) const;
};

template <typename GbxInitConds>
dualview_gbx create_gbxs(const GbxInitConds &gic,
                        const viewd_supers supers)
{

  std::cout << "\n--- create gridboxes ---\n"
            << "initialising\n";
  const dualview_gbx gbxs(initialise_gbxs(gic, supers));

  std::cout << "checking initialisation\n";
  is_gbxinit_complete(gbxs, gic.get_size());
  print_gbxs(gbxs.view_host());

  std::cout << "--- create gridboxes: success ---\n";

  return gbxs;
}

template <typename GbxInitConds>
inline dualview_gbx
initialise_gbxs(const GbxInitConds &gic,
                            const viewd_supers supers)
/* initialise a view of superdrops (on device memory)
using data from an InitData instance for their initial
gbxindex, spatial coordinates and attributes */
{
  // create dualview for gridboxes on device and host memory
  dualview_gbx gbxs("gbxs", gic.get_ngbxs());

  // initialise gridboxes on host
  gbxs.sync_host();
  initialise_gbxs_on_host(gic, supers, gbxs.view_host());
  gbxs.modify_host();

  // update device gridbox view to match host's gridbox view
  gbxs.sync_device();

  return gbxs;
}

template <typename GbxInitConds>
inline void initialise_gbxs_on_host(const GbxInitConds &gic,
                            const viewd_supers supers,
                            const viewh_gbx h_gbxs) 
/* initialise the host view of gridboxes
using some data from a GbxInitConds instance
e.g. for each gridbox's volume */
{
  const size_t ngbxs(h_gbxs.extent(0));
  const GenGridbox gen(gic);

  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    h_gbxs(ii) = gen(ii, supers);
  }
}

#endif // CREATEGBXS_HPP