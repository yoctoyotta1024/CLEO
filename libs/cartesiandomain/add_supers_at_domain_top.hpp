/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: add_supers_at_domain_top.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 18th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * A definition of the Domain Boundary Conditions to use for Cartesian GridBox Maps, Motion of
 * Super-Droplets and MoveSupersInDomain
 */

#ifndef LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
#define LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "gridboxes/sortsupers.hpp"

struct AddSupersAtDomainTop {
  double coord3lim; /**< gridboxes with upper bound > coord3lim get new super-droplets */

  void new_supers_for_gridbox(const Gridbox &gbx, const viewd_supers totsupers) const {
    std::cout << "new supers for gbx: " << gbx.get_gbxindex() << "\n";
  }

  /* (re)sorting supers based on their gbxindexes and then updating the span for each
  gridbox accordingly.
  Kokkos::parallel_for([...]) (on host) is equivalent to:
  for (size_t ii(0); ii < ngbxs; ++ii){[...]}
  when in serial */
  void move_supers_between_gridboxes(const viewd_gbx d_gbxs, const viewd_supers totsupers) const {
    sort_supers(totsupers);

    const size_t ngbxs(d_gbxs.extent(0));
    Kokkos::parallel_for(
        "move_supers_between_gridboxes", TeamPolicy(ngbxs, Kokkos::AUTO()),
        KOKKOS_CLASS_LAMBDA(const TeamMember &team_member) {
          const int ii = team_member.league_rank();

          auto &gbx(d_gbxs(ii));
          gbx.supersingbx.set_refs(team_member);
        });
  }

  /* New super-droplets are added to domain with coord3 >= COORD3LIM [m] */
  explicit AddSupersAtDomainTop(const double COORD3LIM) : coord3lim(COORD3LIM / dlc::COORD0) {}

  void operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                  const viewd_supers totsupers) const {
    const size_t ngbxs(d_gbxs.extent(0));

    bool is_supers_added = false;
    for (size_t ii(0); ii < ngbxs; ++ii) {  // TODO(CB) parallelise?
      auto &gbx(d_gbxs(ii));
      const auto bounds = gbxmaps.coord3bounds(gbx.get_gbxindex());
      if (bounds.second > coord3lim) {
        new_supers_for_gridbox(gbx, totsupers);
        is_supers_added = true;
      }
    }

    if (is_supers_added) {  // resort totsupers view and set gbx references
      move_supers_between_gridboxes(d_gbxs, totsupers);
    }
  }
};

#endif  // LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
