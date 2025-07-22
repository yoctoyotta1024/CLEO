/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cartesian_transport_across_domain.hpp
 * Project: movement
 * Created Date: Monday 24th Febuary 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Movement of a superdroplet throughout a cartesian domain, optionally distributed across
 * more than one MPI process
 */

#ifndef LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_TRANSPORT_ACROSS_DOMAIN_HPP_
#define LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_TRANSPORT_ACROSS_DOMAIN_HPP_

#include <Kokkos_Core.hpp>
#include <array>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

#include "../../cleoconstants.hpp"
#include "../../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "configuration/communicator.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/supersindomain.hpp"
#include "mpi.h"
#include "superdrops/superdrop.hpp"

/*
function to move super-droplets between MPI processes, e.g. for superdroplets
which move to/from gridboxes on different nodes.
*/
template <GridboxMaps GbxMaps>
viewd_supers sendrecv_supers(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                             viewd_supers totsupers);

/*
 * struct satisfying TransportAcrossDomain concept for transporting superdroplets around a
 * cartesian domain, optionally with MPI communicaiton of superdroplets between nodes
 */
struct CartesianTransportAcrossDomain {
  /* (re)sorting supers based on their gbxindexes as step to 'move' superdroplets across the domain.
  May also include MPI communication with moves superdroplets away from/into a node's domain
  */
  SupersInDomain operator()(const CartesianMaps &gbxmaps, const viewd_gbx d_gbxs,
                            SupersInDomain &allsupers) const;
};

/*
function to move super-droplets between MPI processes, e.g. for superdroplets
which move to/from gridboxes on different nodes.
*/
template <GridboxMaps GbxMaps>
viewd_supers sendrecv_supers(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                             viewd_supers totsupers) {
  int comm_size, my_rank;
  comm_size = init_communicator::get_comm_size();
  my_rank = init_communicator::get_comm_rank();

  std::vector<MPI_Request> exchange_requests(comm_size * 6, MPI_REQUEST_NULL);
  std::vector<MPI_Status> exchange_statuses(comm_size * 6);

  std::vector<int> per_process_send_superdrops(comm_size, 0);
  std::vector<int> per_process_recv_superdrops(comm_size, 0);

  std::vector<int> uint64_send_displacements(comm_size, 0);
  std::vector<int> uint64_recv_displacements(comm_size, 0);

  std::vector<int> uint_send_displacements(comm_size, 0);
  std::vector<int> uint_recv_displacements(comm_size, 0);
  std::vector<int> uint_send_counts(comm_size, 0);
  std::vector<int> uint_recv_counts(comm_size, 0);

  std::vector<int> double_send_displacements(comm_size, 0);
  std::vector<int> double_recv_displacements(comm_size, 0);
  std::vector<int> double_send_counts(comm_size, 0);
  std::vector<int> double_recv_counts(comm_size, 0);
  std::vector<std::vector<int>> superdrops_indices_per_process(comm_size);

  size_t total_superdrops_to_send = 0;
  size_t total_superdrops_to_recv = 0;
  size_t local_superdrops = 0;
  size_t superdrop_index = totsupers.extent(0) - 1;
  Superdrop &drop = totsupers(superdrop_index);

  // Go through superdrops from back to front and find how many should be sent and their indices
  const auto ngbxs = d_gbxs.extent(0);
  while (drop.get_sdgbxindex() >= ngbxs) {
    if (drop.get_sdgbxindex() < LIMITVALUES::oob_gbxindex) {
      int target_process = (LIMITVALUES::oob_gbxindex - drop.get_sdgbxindex()) - 1;
      per_process_send_superdrops[target_process]++;
      superdrops_indices_per_process[target_process].push_back(superdrop_index);
      total_superdrops_to_send++;
    }
    drop = totsupers(--superdrop_index);
  }
  local_superdrops = superdrop_index + 1;

  // Share how many superdrops each process will send and receive to/from the others
  MPI_Alltoall(per_process_send_superdrops.data(), 1, MPI_INT, per_process_recv_superdrops.data(),
               1, MPI_INT, MPI_COMM_WORLD);
  total_superdrops_to_recv =
      std::accumulate(per_process_recv_superdrops.begin(), per_process_recv_superdrops.end(), 0);

  assert((local_superdrops + total_superdrops_to_recv <= totsupers.extent(0)) &&
         "must have enough space in supers view to receive superdroplets");
  if (local_superdrops + total_superdrops_to_recv > totsupers.extent(0)) {
    std::cout << "MAXIMUM NUMBER OF LOCAL SUPERDROPLETS EXCEEDED" << std::endl;
    return totsupers;
  }

  // Knowing how many superdroplets will be sent and received, allocate
  // buffers to serialize the data
  std::vector<double> superdrops_double_send_data(total_superdrops_to_send * 5);
  std::vector<double> superdrops_double_recv_data(total_superdrops_to_recv * 5);
  std::vector<unsigned int> superdrops_uint_send_data(total_superdrops_to_send * 2);
  std::vector<unsigned int> superdrops_uint_recv_data(total_superdrops_to_recv * 2);
  std::vector<uint64_t> superdrops_uint64_send_data(total_superdrops_to_send);
  std::vector<uint64_t> superdrops_uint64_recv_data(total_superdrops_to_recv);

  // Calculate the send and receive counts and displacements for each of the target processes
  for (int i = 0; i < comm_size; i++) {
    double_send_counts[i] = per_process_send_superdrops[i] * 5;
    double_recv_counts[i] = per_process_recv_superdrops[i] * 5;
    uint_send_counts[i] = per_process_send_superdrops[i] * 2;
    uint_recv_counts[i] = per_process_recv_superdrops[i] * 2;
    if (i > 0) {
      uint_send_displacements[i] = uint_send_displacements[i - 1] + uint_send_counts[i - 1];
      uint_recv_displacements[i] = uint_recv_displacements[i - 1] + uint_recv_counts[i - 1];

      uint64_send_displacements[i] =
          uint64_send_displacements[i - 1] + per_process_send_superdrops[i - 1];
      uint64_recv_displacements[i] =
          uint64_recv_displacements[i - 1] + per_process_recv_superdrops[i - 1];

      double_send_displacements[i] = double_send_displacements[i - 1] + double_send_counts[i - 1];
      double_recv_displacements[i] = double_recv_displacements[i - 1] + double_recv_counts[i - 1];
    }
  }

  // Serialize the data for all superdroplets into the exchange arrays
  unsigned int send_superdrop_index = 0;
  for (int process_index = 0; process_index < comm_size; process_index++)
    for (int superdrop = 0; superdrop < per_process_send_superdrops[process_index]; superdrop++) {
      superdrop_index = superdrops_indices_per_process[process_index][superdrop];
      totsupers(superdrop_index)
          .serialize_uint_components(superdrops_uint_send_data.begin() + send_superdrop_index * 2);
      totsupers(superdrop_index)
          .serialize_uint64_components(superdrops_uint64_send_data.begin() + send_superdrop_index);
      totsupers(superdrop_index)
          .serialize_double_components(superdrops_double_send_data.begin() +
                                       send_superdrop_index * 5);
      send_superdrop_index++;
    }

  for (int i = 0; i < comm_size; i++) {
    if (i != my_rank) {
      // Checks whether something should be sent to process i
      if (per_process_send_superdrops[i] > 0) {
        MPI_Isend(superdrops_uint_send_data.data() + uint_send_displacements[i],
                  uint_send_counts[i], MPI_UNSIGNED, i, 0, MPI_COMM_WORLD,
                  exchange_requests.data() + i);

        MPI_Isend(superdrops_uint64_send_data.data() + uint64_send_displacements[i],
                  per_process_send_superdrops[i], MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD,
                  exchange_requests.data() + comm_size + i);

        MPI_Isend(superdrops_double_send_data.data() + double_send_displacements[i],
                  double_send_counts[i], MPI_DOUBLE, i, 2, MPI_COMM_WORLD,
                  exchange_requests.data() + comm_size * 2 + i);
      }

      // Checks whether something should be received from process i
      if (per_process_recv_superdrops[i] > 0) {
        MPI_Irecv(superdrops_uint_recv_data.data() + uint_recv_displacements[i],
                  uint_recv_counts[i], MPI_UNSIGNED, i, 0, MPI_COMM_WORLD,
                  exchange_requests.data() + (comm_size * 3) + i);

        MPI_Irecv(superdrops_uint64_recv_data.data() + uint64_recv_displacements[i],
                  per_process_recv_superdrops[i], MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD,
                  exchange_requests.data() + (comm_size * 4) + i);

        MPI_Irecv(superdrops_double_recv_data.data() + double_recv_displacements[i],
                  double_recv_counts[i], MPI_DOUBLE, i, 2, MPI_COMM_WORLD,
                  exchange_requests.data() + (comm_size * 5) + i);
      }
    }
  }

  MPI_Waitall(comm_size * 6, exchange_requests.data(), exchange_statuses.data());

  for (unsigned int i = local_superdrops; i < local_superdrops + total_superdrops_to_recv; i++) {
    int data_offset = i - local_superdrops;
    totsupers(i).deserialize_components(superdrops_uint_recv_data.begin() + data_offset * 2,
                                        superdrops_uint64_recv_data.begin() + data_offset,
                                        superdrops_double_recv_data.begin() + data_offset * 5);

    // Get the local gridbox index which contains the superdroplet
    auto drop_coords = std::array<double, 3>{totsupers(i).get_coord3(), totsupers(i).get_coord1(),
                                             totsupers(i).get_coord2()};
    const auto b4 = std::array<double, 3>{drop_coords[0], drop_coords[1], drop_coords[2]};
    const auto gbxindex =
        (unsigned int)gbxmaps.get_domain_decomposition().get_local_bounding_gridbox_index(
            drop_coords);  // TODO(ALL): access through gbxmaps (note error in conversions?)

    // Since the coordinates have already been corrected in the sending
    // process here just the gridbox index update is necessary
    assert((drop_coords[0] == b4[0]) && (drop_coords[1] == b4[1]) && (drop_coords[2] == b4[2]) &&
           "drop coordinates should have already been corrected and so shoudn't have changed here");
    totsupers(i).set_sdgbxindex(gbxindex);
  }

  // Reset all remaining non-used superdroplet spots
  for (unsigned int i = local_superdrops + total_superdrops_to_recv; i < totsupers.extent(0); i++)
    totsupers(i).set_sdgbxindex(LIMITVALUES::oob_gbxindex);

  return totsupers;
}

#endif  // LIBS_CARTESIANDOMAIN_MOVEMENT_CARTESIAN_TRANSPORT_ACROSS_DOMAIN_HPP_
