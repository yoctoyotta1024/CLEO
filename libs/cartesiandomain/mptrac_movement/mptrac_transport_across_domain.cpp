/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: mptrac_transport_across_domain.cpp
 * Project: mptrac_movement
 * Created Date: Monday 24th Febuary 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors: Jan Clemens (JC)
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
extern "C" {
  #include "mptrac.h"
}

/* (re)sorting supers based on their gbxindexes as step to 'move' superdroplets across the domain.
May also include MPI communication with moves superdroplets away from/into a node's domain
*/
SupersInDomain MPTRACTransportAcrossDomain::operator()(const CartesianMaps &gbxmaps,
                                                          const viewd_gbx d_gbxs,
                                                          SupersInDomain &allsupers) const
{
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

  std::cout << "FORCE TERMINATE DEBUGGING STOP" << std::flush;
  throw std::runtime_error("DEBUGGING STOP");

  return allsupers;
}

/*
function to move super-droplets between MPI processes using MPTRAC library,
e.g. for superdroplets which move to/from gridboxes on different nodes.
*/
viewd_supers MPTRACTransportAcrossDomain::sendrecv_supers(const CartesianMaps &gbxmaps, const viewd_gbx d_gbxs,
viewd_supers totsupers) const
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
 
  if (comm_rank == 0){
    std::cout << "----------> start MPTRAC <----------\n";
  }

  /* define destination ranks for current MPI rank */
  constexpr int ndestinations = 26;
  auto destinations = std::array<int, ndestinations>{};
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

  /* set pointers in particle_ptr_t type for each superdroplet to superdroplet
  attributes and set value for destination rank of each particle */
  constexpr int nquantities = 8; // must match number of superdroplet data (sdgbxindex -> sdId.value)
  const auto ngbxs = gbxmaps.get_local_ngridboxes(); 
  const auto nparticles = totsupers.extent(0);
  auto particles_ptr = std::vector<particle_ptr_t>{};
  auto target_ranks = std::vector<int>{};
  int is_first_rank = 0;
  size_t is_first_sdId = 0;
  for (size_t ip = 0; ip < nparticles; ++ip){
    auto &drop = totsupers(ip);
    const auto p_ptr = drop._make_mptrac_compatible_ptr<particle_ptr_t>();
    particles_ptr.push_back(p_ptr);

    const auto sdgbxindex = drop.get_sdgbxindex();
    int rank_dest = int{comm_rank}; // default is indicator to MPTRAC that particles stay on current rank
    if (sdgbxindex == LIMITVALUES::oob_gbxindex)
    {
      rank_dest = -1; // particles stay on current rank and can be overwritten
    } else if (!(sdgbxindex < ngbxs))
    {
      // check ```sdgbxindex < ngbxs``` assumes 0 <= local_gbxindexes < ngbxs
      // (!) calculation of drop's rank must match calculation according to domain_decomposition
      rank_dest = (LIMITVALUES::oob_gbxindex - 1) - sdgbxindex;

      if (comm_rank == 0 && is_first_rank == 0){
        std::cout << "first superdrop to move: ip=" << ip << "\n";
        std::cout << "first superdrop to move: rank=" << rank_dest << "\n";
        std::cout << "first superdrop to move: sdId=" << drop.sdId.get_value() << "\n";
        std::cout << "first superdrop to move: q[0]=" << sdgbxindex << "\n";
        std::cout << "first superdrop to move: q[1]=" << *particles_ptr.at(ip).q[1] << "\n";
        std::cout << "first superdrop to move: q[2]=" << *particles_ptr.at(ip).q[2] << "\n";
        std::cout << "first superdrop to move: q[3]=" << *particles_ptr.at(ip).q[3] << "\n";
        std::cout << "first superdrop to move: q[4]=" << *(uint64_t*)(particles_ptr.at(ip).q[4]) << "\n";
        std::cout << "first superdrop to move: q[5]=" << *particles_ptr.at(ip).q[5] << "\n";
        std::cout << "first superdrop to move: q[6]=" << *particles_ptr.at(ip).q[6] << "\n";
        std::cout << "first superdrop to move: q[7]=" << *(uint64_t*)(particles_ptr.at(ip).q[7]) << "\n";
        is_first_rank = rank_dest;
        is_first_sdId = drop.sdId.get_value();
      }

    }
    target_ranks.push_back(rank_dest);
  }
  std::cout << "bcast before: " << is_first_rank << ", " << is_first_sdId << "\n";
  MPI_Bcast(&is_first_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&is_first_sdId, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  std::cout << "bcast after: " << is_first_rank << ", " << is_first_sdId << "\n";

  /* define MPI_Datatype for particle */
  MPI_Datatype MPI_Particle;
  register_MPI_type_particle(&MPI_Particle, nquantities);

  /* call MPTRAC MPI communication of particles */
  auto first_particles_ptr = &particles_ptr.front();
  auto q_sizes = std::array<size_t, nquantities>{4, 8, 8, 8, 8, 8, 8, 8}; // bytes of each item in q
  dd_communicate_particles_cleo(
    first_particles_ptr,
    nparticles,
    MPI_Particle,
    destinations.data(),
    ndestinations,
    target_ranks.data(),
    q_sizes.data()
  );

  if (comm_rank == is_first_rank) {
    for (size_t ip = 0; ip < nparticles; ++ip){
      auto &drop = totsupers(ip);
      if (drop.sdId.get_value() == is_first_sdId){
          std::cout << "first superdrop that moved: dest_ip=" << ip << "\n";
          std::cout << "first superdrop that moved: dest_rank=" << is_first_rank << "\n";
          std::cout << "first superdrop that moved: sdId=" << is_first_sdId << "\n";
          std::cout << "first superdrop that moved: q[0]=" << drop.get_sdgbxindex() << "\n";
          std::cout << "first superdrop that moved: q[1]=" << drop.get_coord3() << "\n";
          std::cout << "first superdrop that moved: q[2]=" << drop.get_coord1() << "\n";
          std::cout << "first superdrop that moved: q[3]=" << drop.get_coord2() << "\n";
          std::cout << "first superdrop that moved: q[4]=" << drop.get_xi() << "\n";
          std::cout << "first superdrop that moved: q[5]=" << drop.get_radius() << "\n";
          std::cout << "first superdrop that moved: q[6]=" << drop.get_msol() << "\n";
          std::cout << "first superdrop that moved: q[7]=" << drop.sdId.get_value() << "\n";
      }
    }
  }

  /* add correction to sdgbxindexes of superdroplets that were sent/recieved */
  for (size_t ip = 0; ip < nparticles; ++ip){
    auto &drop = totsupers(ip);
    const auto sdgbxindex = drop.get_sdgbxindex();
    if (!(sdgbxindex < ngbxs) && sdgbxindex != LIMITVALUES::oob_gbxindex)
    {
      // check ```sdgbxindex < ngbxs``` assumes 0 <= local_gbxindexes < ngbxs
      // if superdroplet sdgbxindex is not already local or out of bounds
      // the superdroplet has been sent or received and therefore sdgbxindex needs correcting

      // (!) calculation of drop's rank must match calculation according to domain_decomposition
      auto drop_rank = (LIMITVALUES::oob_gbxindex - 1) - sdgbxindex; // rank where drop should exist

      if (drop_rank != (unsigned int)comm_rank)
      {
        // drop has been sent so needs to be removed
        drop.set_sdgbxindex(LIMITVALUES::oob_gbxindex);
      } else
      {
        // drop has been received so sdgbxindex needs to be corrected to a local one
        auto drop_coords = std::array<double, 3>{drop.get_coord3(), drop.get_coord1(),
          drop.get_coord2()};
        const auto b4 = std::array<double, 3>{drop_coords[0], drop_coords[1], drop_coords[2]};
        const auto gbxindex = gbxmaps.get_domain_decomposition().get_local_bounding_gridbox(drop_coords);
        drop.set_sdgbxindex(gbxindex);
        assert((drop_coords[0] == b4[0]) && (drop_coords[1] == b4[1]) && (drop_coords[2] == b4[2]) &&
            "drop coordinates should have already been corrected and so shoudn't have changed here");
      }
    }
  }

  if (comm_rank == 0){
    std::cout << "----------> finish MPTRAC <----------\n";
  }

  return totsupers;
}

/* define MPI_Datatype for MPTRAC superdroplet particle */
void MPTRACTransportAcrossDomain::register_MPI_type_particle(MPI_Datatype *MPI_Particle,
  const int nquantities) const
{
  constexpr int nblocks = 1;
  MPI_Datatype types[nblocks] = { MPI_DOUBLE }; // (casting all superdroplet quantities to double first)
  int blocklengths[nblocks] = { nquantities }; // number of superdroplet quantities to send in 'q'
  MPI_Aint displacements[nblocks] = { offsetof (particle_quant_t, q) };
  MPI_Type_create_struct (nblocks, blocklengths, displacements, types, MPI_Particle);
  MPI_Type_commit(MPI_Particle);
}