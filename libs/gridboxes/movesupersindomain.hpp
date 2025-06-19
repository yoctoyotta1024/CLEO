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
 * Last Modified: Wednesday 28th May 2025
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

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "configuration/communicator.hpp"
#include "gridboxes/boundary_conditions.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/supersindomain.hpp"
#include "gridboxes/transport_across_domain.hpp"
#include "mpi.h"
#include "superdrops/motion.hpp"
#include "superdrops/sdmmonitor.hpp"

namespace KCS = KokkosCleoSettings;

/*
MoveSupersInGridboxesFunctor struct encapsulates superdroplet motion so that parallel loop
functors only captures motion and not other members of MoveSupersInDomain coincidentally
(which e.g. may not be GPU compatible).
*/
template <GridboxMaps GbxMaps, Motion<GbxMaps> M>
struct MoveSupersInGridboxesFunctor {
  const M sdmotion;
  const GbxMaps gbxmaps;
  const viewd_gbx d_gbxs;
  const subviewd_supers domainsupers;

  /*
   * enact steps (1) and (2) movement of superdroplets for 1 gridbox:
   * (1) update their spatial coords according to type of sdmotion. (device)
   * (2) update their sdgbxindex accordingly (device).
   *
   * Kokkos::parallel_for([...]) is equivalent to:
   * for (size_t kk(0); kk < supers.extent(0); ++kk) {[...]}
   * when in serial
   */
  KOKKOS_INLINE_FUNCTION
  void move_supers_in_gbx(const TeamMember &team_member, const unsigned int gbxindex,
                          const State &state, const subviewd_supers supers) const {
    const size_t nsupers(supers.extent(0));
    auto &_sdmotion = this->sdmotion;
    auto &_gbxmaps = this->gbxmaps;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, nsupers),
                         [=, &_sdmotion, &_gbxmaps](const size_t kk) {
                           /* step (1) */
                           _sdmotion.superdrop_coords(gbxindex, _gbxmaps, state, supers(kk));

                           /* step (2) */
                           _sdmotion.superdrop_gbx(gbxindex, _gbxmaps, supers(kk));
                         });
  }

  /*
   * operator for functor with parallel (TeamPolicy) loop over gridboxes in
   * d_gbxs view in order to call move_supers_in_gbx
   */
  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamMember &team_member) const {
    const auto ii = team_member.league_rank();
    auto &gbx = d_gbxs(ii);
    move_supers_in_gbx(team_member, gbx.get_gbxindex(), gbx.state, gbx.supersingbx(domainsupers));
  }
};

/*
struct for functionality to move superdroplets throughtout the domain by updating their
spatial coordinates (according to some type of Motion) and then updating their gridbox indexes
concordantly and moving them between gridboxes afterwards (MPI communication and sorting) and
finally applying any additional bounday conditions.
*/
template <GridboxMaps GbxMaps, Motion<GbxMaps> M, TransportAcrossDomain<GbxMaps> T,
          BoundaryConditions<GbxMaps> BCs>
class MoveSupersInDomain {
 public:
  /* (expensive!) test if superdrops' gbxindex doesn't match gridbox's gbxindex,
  raise error is assertion fails */
  void check_sdgbxindex_during_motion(const viewd_constgbx d_gbxs,
                                      const viewd_constsupers totsupers) const {
    const auto ngbxs = d_gbxs.extent(0);
    Kokkos::parallel_for(
        "check_sdgbxindex_during_motion", TeamPolicy(ngbxs, KCS::team_size),
        KOKKOS_LAMBDA(const TeamMember &team_member) {
          const auto ii = team_member.league_rank();
          assert(d_gbxs(ii).supersingbx.iscorrect(team_member, totsupers) &&
                 "incorrect references to superdrops in gridbox during motion");
        });
  }

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
        "set_gridboxes_refs", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
        KOKKOS_LAMBDA(const size_t ii) { d_gbxs(ii).supersingbx.set_refs(domainsupers); });
  }

  /*
   * call operator of MoveSupersInGridboxesFunctor struct to enact steps (1) and (2) of superdroplet
   * motion:
   * (1) update superdroplets' spatial coords according to type of sdmotion. (device)
   * (2) update superdroplets' sdgbxindex accordingly (device).
   *
   * Kokkos::parallel_for([...]) is equivalent to:
   * for (size_t ii(0); ii < ngbxs; ++ii) {[...]}
   * when in serial
   */
  void move_supers_in_gridboxes(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                                const subviewd_supers domainsupers) const {
    Kokkos::Profiling::ScopedRegion region("sdm_movement_move_in_gridboxes");

    const size_t ngbxs(d_gbxs.extent(0));
    const auto functor =
        MoveSupersInGridboxesFunctor<GbxMaps, M>{sdmotion, gbxmaps, d_gbxs, domainsupers};
    Kokkos::parallel_for("move_supers_in_gridboxes", TeamPolicy(ngbxs, KCS::team_size), functor);
  }

 private:
  M sdmotion;
  T transport_across_domain;
  BCs boundary_conditions;

  /*
   * step (3) move superdroplets between gridboxes
   *
   * (re)sorting supers, based on their gbxindexes and then updating the refs for each gridbox
   * accordingly. May also include MPI communication with moves superdroplets away from/into a
   * node's domain and/or check that superdroplets end up in correct gridboxes.
   *
   * optionally can incude test to check superdroplets end up in correct gridbox (for debugging)
   */
  SupersInDomain move_supers_between_gridboxes(const GbxMaps &gbxmaps, const viewd_gbx d_gbxs,
                                               SupersInDomain &allsupers) const {
    Kokkos::Profiling::ScopedRegion region("sdm_movement_between_gridboxes");

    allsupers = transport_across_domain(gbxmaps, d_gbxs, allsupers);

    set_gridboxes_refs(d_gbxs, allsupers.domain_supers());

    /* optional (expensive!) test if superdrops' gbxindex doesn't match gridbox's gbxindex */
    // check_sdgbxindex_during_motion(d_gbxs, allsupers.get_totsupers_readonly());

    return allsupers;
  }

  /* enact movement of superdroplets throughout domain in three stages:
  (1) update their spatial coords according to type of sdmotion. (device)
  (2) update their sdgbxindex accordingly (device)
  (3) move superdroplets between gridboxes (host)
  (4) apply domain boundary conditions (host and/or device)
  */
  SupersInDomain move_superdrops_in_domain(const unsigned int t_sdm, const GbxMaps &gbxmaps,
                                           viewd_gbx d_gbxs, SupersInDomain &allsupers) const {
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    /* steps (1 - 2) */
    move_supers_in_gridboxes(gbxmaps, d_gbxs, allsupers.domain_supers());

    /* step (3) */
    allsupers = move_supers_between_gridboxes(gbxmaps, d_gbxs, allsupers);

    /* step (4) */
    Kokkos::Profiling::pushRegion("sdm_movement_boundary_conditions");
    allsupers = boundary_conditions.apply(gbxmaps, d_gbxs, allsupers);
    Kokkos::Profiling::popRegion();

    return allsupers;
  }

 public:
  MoveSupersInDomain(const M i_sdmotion, const T i_transport_across_domain,
                     const BCs i_boundary_conditions)
      : sdmotion(i_sdmotion),
        transport_across_domain(i_transport_across_domain),
        boundary_conditions(i_boundary_conditions) {}

  /* extra constructor useful to help when compiler cannot deduce type of GbxMaps */
  MoveSupersInDomain(const GbxMaps &gbxmaps, const M i_sdmotion, const T i_transport_across_domain,
                     const BCs i_boundary_conditions)
      : MoveSupersInDomain(i_sdmotion, i_transport_across_domain, i_boundary_conditions) {}

  /* returns time when superdroplet motion is next due to occur given current time, t_sdm */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const { return sdmotion.next_step(t_sdm); }

  /*
   * if current time, t_sdm, is time when superdrop motion should occur, enact movement of
   * superdroplets throughout domain.
   *
   * @param allsupers Struct to handle all superdrops (both in and out of bounds of domain).
   *
   */
  SupersInDomain run_step(const unsigned int t_sdm, const GbxMaps &gbxmaps, viewd_gbx d_gbxs,
                          SupersInDomain &allsupers, const SDMMonitor auto mo) const {
    if (sdmotion.on_step(t_sdm)) {
      allsupers = move_superdrops_in_domain(t_sdm, gbxmaps, d_gbxs, allsupers);
      mo.monitor_motion(d_gbxs, allsupers.domain_supers_readonly());
    }

    return allsupers;
  }
};

#endif  // LIBS_GRIDBOXES_MOVESUPERSINDOMAIN_HPP_
