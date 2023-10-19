/*
 * ----- CLEO -----
 * File: creategbxs.hpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 19th October 2023
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
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/state.hpp"

class CreateGbxs
{
private:
  class GenGridbox
  {
  private:
    std::unique_ptr<Gridbox::Gbxindex::Gen> GbxindexGen; // pointer to gridbox index generator
    std::vector<double> volumes;

    inline State state_at(const unsigned int ii) const;
  
  public:
    template <typename FetchInitData>
    inline GenGridbox(const FetchInitData &fid);

    inline Gridbox operator()(const unsigned int ii,
                              const viewd_constsupers supers) const;
  };

  template <typename FetchInitData>
  inline dualview_gbx initialise_gbxs(const FetchInitData &fid,
                                      const viewd_constsupers supers) const;
  /* initialise a dualview of gridboxes (on host and device 
  memory) using data from a FetchInitData instance to initialise
  the host view and then syncing the view to the device */

  template <typename FetchInitData>
  inline void initialise_gbxs_on_host(const FetchInitData &fid,
                                      const viewd_constsupers supers,
                                      const viewh_gbx h_gbxs) const;
  /* initialise the host view of gridboxes
  using some data from a FetchInitData instance
  e.g. for each gridbox's volume */

  void ensure_initialisation_complete(dualview_gbx gbxs,
                                      const size_t size) const;

  void print_gbxs(const viewh_constgbx gbxs) const;
  /* print gridboxes information */

public:
  template <typename FetchInitData>
  dualview_gbx operator()(const FetchInitData fid,
                          const viewd_constsupers supers) const
  {

    std::cout << "\n--- create gridboxes ---\n"
              << "initialising\n";
    const dualview_gbx gbxs(initialise_gbxs(fid, supers));
 
    std::cout << "checking initialisation\n";
    ensure_initialisation_complete(gbxs, fid.get_size());
    print_gbxs(gbxs.view_host());

    std::cout << "--- create gridboxes: success ---\n";

    return gbxs;
  }
};

template <typename FetchInitData>
inline dualview_gbx
CreateGbxs::initialise_gbxs(const FetchInitData &fid,
                            const viewd_constsupers supers) const
/* initialise a view of superdrops (on device memory)
using data from an InitData instance for their initial
gbxindex, spatial coordinates and attributes */
{
  // create dualview for gridboxes on device and host memory
  dualview_gbx gbxs("gbxs", fid.get_ngbxs());

  // initialise gridboxes on host
  gbxs.sync_host();
  initialise_gbxs_on_host(fid, supers, gbxs.view_host());
  gbxs.modify_host();

  // update device gridbox view to match host's gridbox view
  gbxs.sync_device();

  return gbxs;
}

template <typename FetchInitData>
inline void CreateGbxs::
    initialise_gbxs_on_host(const FetchInitData &fid,
                            const viewd_constsupers supers,
                            const viewh_gbx h_gbxs) const
/* initialise the host view of gridboxes
using some data from a FetchInitData instance
e.g. for each gridbox's volume */
{
  const size_t ngbxs(h_gbxs.extent(0));
  const GenGridbox gen_gridbox(fid);

  for (size_t ii(0); ii < ngbxs; ++ii)
  {
    h_gbxs(ii) = gen_gridbox(ii, supers);
  }
}

template <typename FetchInitData>
inline CreateGbxs::GenGridbox::
    GenGridbox(const FetchInitData &fid)
    : GbxindexGen(std::make_unique<Gridbox::Gbxindex::Gen>()),
      volumes(fid.volume()) {}

inline Gridbox
CreateGbxs::GenGridbox::operator()(const unsigned int ii,
                                   const viewd_constsupers supers) const
{
  const auto gbxindex(GbxindexGen->next());
  const State state(state_at(ii)); 
  
  return Gridbox(gbxindex, state, supers);
}

inline State 
CreateGbxs::GenGridbox::state_at(const unsigned int ii) const
{
  double volume(volumes.at(ii));
  double press(0.0);                   
  double temp(0.0);                    
  double qvap(0.0);                    
  double qcond(0.0);                   
  Kokkos::pair<double, double> wvel{0.0,0.0}; 
  Kokkos::pair<double, double> uvel{0.0,0.0};
  Kokkos::pair<double, double> vvel{0.0,0.0};

  return State(volume,
               press, temp, qvap, qcond,
               wvel, uvel, vvel);
}

#endif // CREATEGBXS_HPP