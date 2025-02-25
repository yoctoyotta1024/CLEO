/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: mptrac_transport_across_domain.cpp
 * Project: mptrac_movement
 * Created Date: Monday 24th Febuary 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 24th Febuary 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Implementation file for movement of a superdroplet throughout a cartesian domain
 * using the MPTRAC library
 */

#include "./mptrac_transport_across_domain.hpp"

/* (re)sorting supers based on their gbxindexes as step to 'move' superdroplets across the domain.
May also include MPI communication with moves superdroplets away from/into a node's domain
*/
SupersInDomain MPTRACTransportAcrossDomain::operator()(const CartesianMaps &gbxmaps,
                                                          const viewd_gbx d_gbxs,
                                                          SupersInDomain &allsupers) const {
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // TODO(ALL): remove guard once domain decomposition is GPU compatible
  if (comm_size > 1) {
    // TODO(ALL): combine two sorts into one(?)
    auto totsupers = allsupers.sort_totsupers_without_set(d_gbxs);
    totsupers = sendrecv_supers(gbxmaps, d_gbxs, totsupers); 
    allsupers.sort_and_set_totsupers(totsupers, d_gbxs);
  } else {
    allsupers.sort_totsupers(d_gbxs);
  }

  return allsupers;
}

/*
function to move super-droplets between MPI processes using MPTRAC library,
e.g. for superdroplets which move to/from gridboxes on different nodes.
*/
 viewd_supers MPTRACTransportAcrossDomain::sendrecv_supers(const CartesianMaps &gbxmaps, const viewd_gbx d_gbxs,
  viewd_supers totsupers) const {

    /* define destination ranks for current MPI rank */
    auto destinations = std::array<int, 26>{};
    const auto neighboring_processes = gbxmaps.get_domain_decomposition().get_neighboring_processes();
    auto dd = 0;
    for (auto k : {-1, 0, 1}){
      for (auto i : {-1, 0, 1}){
        for (auto j : {-1, 0, 1}) {
          if (i != 0 || j != 0 || k != 0) {
            destinations[dd] = neighboring_processes.at({k, i, j});
            ++dd;
          }
        }
      }
    }

    /* set pointers in particle_ptr_t type for each superdroplet to superdroplet attributes */
   const auto nparticles = totsupers.extent(0);
    auto particles_ptr = std::vector<std::shared_ptr<particle_ptr_t>>{};
    for (size_t ip = 0; ip < nparticles; ++ip){
      auto &drop = totsupers(ip);
      const auto p_ptr = drop._make_mptrac_compatible_ptr<particle_ptr_t>();
      particles_ptr.push_back(p_ptr);
    }

    return totsupers;
  }