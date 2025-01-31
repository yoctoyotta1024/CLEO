/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: movesupersindomain.hpp
 * Project: gridboxes
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality related to moving superdroplets
 * (both updating their spatial coordinates and
 * moving them between gridboxes)
 */

#ifndef LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_
#define LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <numeric>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/supersindomain.hpp"
#include "mpi.h"
#include "superdrops/motion.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "superdrops/superdrop.hpp"

/*
function to move super-droplets between MPI processes, e.g. for superdroplets
which move to/from gridboxes on different nodes.
*/
template <GridboxMaps GbxMaps>
viewd_supers sendrecv_supers(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                             viewd_supers totsupers);

/*
struct for functionality to move superdroplets throughtout
the domain by updating their spatial coordinates (according to
some type of Motion) and then moving them between gridboxes
after updating their gridbox indexes concordantly
*/
template <GridboxMaps GbxMaps, Motion<GbxMaps> M, typename BoundaryConditions>
struct MoveSupersInDomain {
  /*
  EnactMotion struct encapsulates motion so that parallel loops with KOKKOS_CLASS_LAMBDA
  (ie. [=] on CPUs) functors only captures motion and not other members of MoveSupersInDomain
  coincidentally (which may not be GPU compatible).
  */
  struct EnactMotion {
    M motion;

    /* enact steps (1) and (2) movement of superdroplets for 1 gridbox:
    (1) update their spatial coords according to type of motion. (device)
    (2) update their sdgbxindex accordingly (device).
    Kokkos::parallel_for([...]) is equivalent to:
    for (size_t kk(0); kk < supers.extent(0); ++kk) {[...]}
    when in serial */
    KOKKOS_INLINE_FUNCTION
    void move_supers_in_gbx(const TeamMember &team_member, const unsigned int gbxindex,
                            const GbxMaps &gbxmaps, const State &state,
                            const subviewd_supers supers) const {
      const size_t nsupers(supers.extent(0));
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, nsupers),
                           [&, this](const size_t kk) {
                             /* step (1) */
                             motion.superdrop_coords(gbxindex, gbxmaps, state, supers(kk));

                             /* step (2) */
                             motion.superdrop_gbx(gbxindex, gbxmaps, supers(kk));
                           });
    }

    /* enact steps (1) and (2) movement of superdroplets
    throughout domain (i.e. for all gridboxes):
    (1) update their spatial coords according to type of motion. (device)
    (2) update their sdgbxindex accordingly (device).
    Kokkos::parallel_for([...]) is equivalent to:
    for (size_t ii(0); ii < ngbxs; ++ii) {[...]}
    when in serial */
    void move_supers_in_gridboxes(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                                  const subviewd_supers domainsupers) const {
      Kokkos::Profiling::ScopedRegion region("sdm_movement_move_in_gridboxes");

      const size_t ngbxs(d_gbxs.extent(0));
      Kokkos::parallel_for(
          "move_supers_in_gridboxes", TeamPolicy(ngbxs, Kokkos::AUTO()),
          KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
            const auto ii = team_member.league_rank();

            auto &gbx = d_gbxs(ii);
            move_supers_in_gbx(team_member, gbx.get_gbxindex(), gbxmaps, gbx.state,
                               gbx.supersingbx(domainsupers));
          });
    }
  } enactmotion;

  /*
  Updates the refs for each gridbox given domainsupers containing all the superdroplets within
  the domain (on one node).
  Kokkos::parallel_for([...]) (on host) is equivalent to:
  for (size_t ii(0); ii < ngbxs; ++ii){[...]}
  when in serial.
  */
  void set_gridboxes_refs(const viewd_gbx d_gbxs, const subviewd_constsupers domainsupers) const {
    const auto ngbxs = d_gbxs.extent(0);
    Kokkos::parallel_for(
        "set_gridboxes_refs", TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_LAMBDA(const TeamMember &team_member) {
          const auto ii = team_member.league_rank();
          d_gbxs(ii).supersingbx.set_refs(team_member, domainsupers);
        });
  }

  /* (expensive!) test if superdrops' gbxindex doesn't match gridbox's gbxindex,
  raise error is assertion fails */
  void check_sdgbxindex_during_motion(const viewd_constgbx d_gbxs,
                                      const viewd_constsupers totsupers) const {
    const auto ngbxs = d_gbxs.extent(0);
    Kokkos::parallel_for(
        "check_sdgbxindex_during_motion", TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_LAMBDA(const TeamMember &team_member) {
          const auto ii = team_member.league_rank();
          assert(d_gbxs(ii).supersingbx.iscorrect(team_member, totsupers) &&
                 "incorrect references to superdrops in gridbox during motion");
        });
  }

  MoveSupersInDomain(const M mtn, const BoundaryConditions boundary_conditions)
      : enactmotion({mtn}), apply_domain_boundary_conditions(boundary_conditions) {}

  /* extra constructor useful to help when compiler cannot deduce type of GBxMaps */
  MoveSupersInDomain(const GbxMaps &gbxmaps, const M mtn,
                     const BoundaryConditions boundary_conditions)
      : MoveSupersInDomain(mtn, boundary_conditions) {}

  /* returns time when superdroplet motion is
  next due to occur given current time, t_sdm */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const {
    return enactmotion.motion.next_step(t_sdm);
  }

  /*
   * if current time, t_sdm, is time when superdrop motion should occur, enact movement of
   * superdroplets throughout domain.
   *
   * @param allsupers Struct to handle all superdrops (both in and out of bounds of domain).
   *
   */
  SupersInDomain run_step(const unsigned int t_sdm, const GbxMaps &gbxmaps, viewd_gbx d_gbxs,
                          SupersInDomain &allsupers, const SDMMonitor auto mo) const {
    if (enactmotion.motion.on_step(t_sdm)) {
      allsupers = move_superdrops_in_domain(t_sdm, gbxmaps, d_gbxs, allsupers);
      mo.monitor_motion(d_gbxs, allsupers.domain_supers_readonly());
    }

    return allsupers;
  }

 private:
  BoundaryConditions apply_domain_boundary_conditions;

  /* (re)sorting supers, based on their gbxindexes and then updating the refs for each gridbox
  accordingly. May also include MPI communication with moves superdroplets away from/into a node's
  domain and/or check that superdroplets end up in correct gridboxes.
  */
  SupersInDomain move_supers_between_gridboxes(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                                               SupersInDomain &allsupers) const {
    Kokkos::Profiling::ScopedRegion region("sdm_movement_between_gridboxes");

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

    set_gridboxes_refs(d_gbxs, allsupers.domain_supers());

    /* optional (expensive!) test if superdrops' gbxindex doesn't match gridbox's gbxindex */
    // check_sdgbxindex_during_motion(d_gbxs, allsupers.get_totsupers_readonly());

    return allsupers;
  }

  /* enact movement of superdroplets throughout domain in three stages:
  (1) update their spatial coords according to type of motion. (device)
  (2) update their sdgbxindex accordingly (device)
  (3) move superdroplets between gridboxes (host)
  (4) (optional) apply domain boundary conditions (host and device)
  // TODO(all) use tasking to convert all 3 team policy
  // loops from first two function calls into 1 loop?
  */
  SupersInDomain move_superdrops_in_domain(const unsigned int t_sdm, const GbxMaps &gbxmaps,
                                           viewd_gbx d_gbxs, SupersInDomain &allsupers) const {
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* steps (1 - 2) */
    enactmotion.move_supers_in_gridboxes(gbxmaps, d_gbxs, allsupers.domain_supers());

    /* step (3) */
    allsupers = move_supers_between_gridboxes(gbxmaps, d_gbxs, allsupers);

    /* step (4) */
    Kokkos::Profiling::pushRegion("sdm_movement_boundary_conditions");
    allsupers = apply_domain_boundary_conditions(gbxmaps, d_gbxs, allsupers);
    Kokkos::Profiling::popRegion();

    return allsupers;
  }
};

/*
function to move super-droplets between MPI processes, e.g. for superdroplets
which move to/from gridboxes on different nodes.
*/
template <GridboxMaps GbxMaps>
viewd_supers sendrecv_supers(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                             viewd_supers totsupers) {
  int comm_size, my_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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
        (unsigned int)gbxmaps.get_domain_decomposition().get_local_bounding_gridbox(
            drop_coords);  // TODO(ALL): access through gbxmaps (note error in conversions?)

    // Since the coordinates have already been corrected in the sending
    // process here just the gridbox index update is necessary
    assert((drop_coords[0] == b4[0]) && (drop_coords[1] == b4[1]) && (drop_coords[2] == b4[2]) &&
           "drop coordinates should have already been corrected and so shoudn't have changed here");
    totsupers(i).set_sdgbxindex(gbxindex);

    // TODO(ALL): add check_bounds to SD?
  }

  // Reset all remaining non-used superdroplet spots
  for (unsigned int i = local_superdrops + total_superdrops_to_recv; i < totsupers.extent(0); i++)
    totsupers(i).set_sdgbxindex(LIMITVALUES::oob_gbxindex);

  return totsupers;
}

#endif  // LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_
