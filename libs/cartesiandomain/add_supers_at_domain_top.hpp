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
#include "cartesiandomain/domainboundaries.hpp"
#include "gridboxes/sortsupers.hpp"
#include "initialise/optional_config_params.hpp"
#include "superdrops/superdrop.hpp"

struct AddSupersAtDomainTop {
  double coord3lim; /**< gridboxes with upper bound > coord3lim get new super-droplets */
  size_t nsupers;   /**< number of superdroplets to add to gridboxes above coord3lim */

  void remove_supers_from_gridbox(const Gridbox &gbx) const {
    const auto supers = gbx.supersingbx();
    for (size_t kk(0); kk < supers.extent(0); ++kk) {
      if (supers(kk).get_coord3() >= coord3lim) {
        supers(kk).set_sdgbxindex(outofbounds_gbxindex());  // remove super-droplet from domain
      }
    }
  }

  // Superdrop create_superdrop(const Gridbox &gbx, const viewd_supers totsupers) const {
  //   const auto sdgbxindex = gbx.get_
  //   return Superdrop(sdgbxindex, coords312[0], coords312[1], coords312[2], attrs, sdId);
  // }

  void add_supers_for_gridbox(const CartesianMaps &gbxmaps, const Gridbox &gbx,
                              const viewd_supers totsupers) const {
    // const auto sd = create_superdrop();
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
  explicit AddSupersAtDomainTop(const OptionalConfigParams::AddSupersAtDomainTopParams &config)
      : coord3lim(config.COORD3LIM / dlc::COORD0), nsupers(config.nsupers) {}

  void operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                  const viewd_supers totsupers) const {
    const size_t ngbxs(d_gbxs.extent(0));

    bool is_supers_added = false;
    for (size_t ii(0); ii < ngbxs; ++ii) {  // TODO(CB) parallelise?
      auto &gbx(d_gbxs(ii));
      const auto bounds = gbxmaps.coord3bounds(gbx.get_gbxindex());
      if (bounds.second > coord3lim) {
        remove_supers_from_gridbox(gbx);
        add_supers_for_gridbox(gbxmaps, gbx, totsupers);
        is_supers_added = true;
      }
    }

    if (is_supers_added) {  // resort totsupers view and set gbx references
      move_supers_between_gridboxes(d_gbxs, totsupers);
    }
  }
};

#endif  // LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
