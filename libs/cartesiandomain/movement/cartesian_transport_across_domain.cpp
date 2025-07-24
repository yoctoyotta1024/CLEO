/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesian_transport_across_domain.cpp
 * Project: movement
 * Created Date: Monday 24th Febuary 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----q
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Implementation file for movement of a superdroplet throughout a cartesian domain,
 * optionally distributed across more than one MPI process
 */

#include "./cartesian_transport_across_domain.hpp"

/* (re)sorting supers based on their gbxindexes as step to 'move' superdroplets across the domain.
May also include MPI communication with moves superdroplets away from/into a node's domain
*/
SupersInDomain CartesianTransportAcrossDomain::operator()(const CartesianMaps &gbxmaps,
                                                          const viewd_gbx d_gbxs,
                                                          SupersInDomain &allsupers) const {
  int comm_size;
  comm_size = init_communicator::get_comm_size();

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
