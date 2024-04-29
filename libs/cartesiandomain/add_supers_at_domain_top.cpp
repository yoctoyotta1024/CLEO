/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: add_supers_at_domain_top.cpp
 * Project: cartesiandomain
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 29th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Implementation of AddSupersAtDomainTop definition of the Domain Boundary Conditions to use
 * for Cartesian GridBox Maps, Motion of Super-Droplets and MoveSupersInDomain
 */

#include "cartesiandomain/add_supers_at_domain_top.hpp"

void AddSupersAtDomainTop::operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
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

/* (re)sorting supers based on their gbxindexes and then updating the span for each
gridbox accordingly.
Kokkos::parallel_for([...]) (on host) is equivalent to:
for (size_t ii(0); ii < ngbxs; ++ii){[...]}
when in serial */
void AddSupersAtDomainTop::move_supers_between_gridboxes(const viewd_gbx d_gbxs,
                                                         const viewd_supers totsupers) const {
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

/* set super-droplet sdgbxindex to out of bounds value */
void AddSupersAtDomainTop::remove_supers_from_gridbox(const Gridbox &gbx) const {
  const auto supers = gbx.supersingbx();
  for (size_t kk(0); kk < supers.extent(0); ++kk) {
    if (supers(kk).get_coord3() >= coord3lim) {
      supers(kk).set_sdgbxindex(outofbounds_gbxindex());  // remove super-droplet from domain
    }
  }
}

void AddSupersAtDomainTop::add_supers_for_gridbox(const CartesianMaps &gbxmaps, const Gridbox &gbx,
                                                  const viewd_supers totsupers) const {
  size_t start = gbx.domain_totnsupers();

  if (start + newnsupers > totsupers.extent(0)) {
    const auto err = std::string(
        "UNDEFINED BEHAVIOUR: Number of super-droplets in the domain cannot become larger than the "
        "size of the super-droplets' view");
    throw std::invalid_argument(err);
  }

  for (size_t kk(start); kk < start + newnsupers; ++kk) {
    totsupers(kk) = create_superdrop(gbx);
  }
}

Superdrop AddSupersAtDomainTop::create_superdrop(const Gridbox &gbx) const {
  const auto sdgbxindex = gbx.get_gbxindex();
  const auto coords312 = create_superdrop_coords();
  const auto attrs = create_superdrop_attrs();
  const auto sdId = sdIdGen->next();

  return Superdrop(sdgbxindex, coords312[0], coords312[1], coords312[2], attrs, sdId);
}

std::array<double, 3> AddSupersAtDomainTop::create_superdrop_coords() const {
  const auto coord3 = double{830.0 / dlc::COORD0};  // TODO(CB): WIP
  const auto coord1 = double{10.0 / dlc::COORD0};
  const auto coord2 = double{10.0 / dlc::COORD0};

  return std::array<double, 3>{coord3, coord1, coord2};
}

SuperdropAttrs AddSupersAtDomainTop::create_superdrop_attrs() const {
  const auto radius = double{1e-3 / dlc::R0};  // TODO(CB): WIP
  const auto msol = double{1e-12 / dlc::MASS0};
  const auto xi = size_t{20};
  const auto solute = SoluteProperties{};

  return SuperdropAttrs(solute, xi, radius, msol);
}
