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
 * Last Modified: Saturday 4th May 2024
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

/* (re)sorting supers based on their gbxindexes and then updating the span for each
gridbox accordingly.
Kokkos::parallel_for([...]) (on host) is equivalent to:
for (size_t ii(0); ii < ngbxs; ++ii){[...]}
when in serial */
void move_supers_between_gridboxes_again(const viewd_gbx d_gbxs, const viewd_supers totsupers) {
  sort_supers(totsupers);

  const size_t ngbxs(d_gbxs.extent(0));
  Kokkos::parallel_for(
      "move_supers_between_gridboxes_again", TeamPolicy(ngbxs, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const TeamMember &team_member) {
        const int ii = team_member.league_rank();

        auto &gbx(d_gbxs(ii));
        gbx.supersingbx.set_refs(team_member);
      });
}

/* set super-droplet sdgbxindex to out of bounds value if superdrop coord3 > coord3lim */
KOKKOS_FUNCTION
void remove_superdrops(const TeamMember &team_member, const Gridbox &gbx, const double coord3lim) {
  const auto supers = gbx.supersingbx();
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team_member, supers.extent(0)), [supers, coord3lim](const size_t kk) {
        if (supers(kk).get_coord3() >= coord3lim) {
          supers(kk).set_sdgbxindex(outofbounds_gbxindex());  // remove super-droplet from domain
        }
      });
}

/* set super-droplet sdgbxindex to out of bounds value if superdrop coord3 > coord3lim */
void remove_superdrops_from_gridboxes(const CartesianMaps &gbxmaps, const viewd_gbx d_gbxs,
                                      const double coord3lim) {
  const size_t ngbxs(d_gbxs.extent(0));
  Kokkos::parallel_for(
      "remove_superdrops", TeamPolicy(ngbxs, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const TeamMember &team_member) {
        const int ii = team_member.league_rank();

        const auto ubound = gbxmaps.coord3bounds(d_gbxs(ii).get_gbxindex()).second;
        if (ubound > coord3lim) {
          remove_superdrops(team_member, d_gbxs(ii), coord3lim);
        }
      });
}

/* create 'newnsupers' number of new superdroplets from the create_superdrop function */
void add_superdrops_for_gridbox(const CartesianMaps &gbxmaps, const Gridbox &gbx,
                                const viewd_supers totsupers) {
  const auto gbxindex = gbx.get_gbxindex();
  const size_t start = gbx.domain_totnsupers();

  if (start + newnsupers > totsupers.extent(0)) {
    const auto err = std::string(
        "UNDEFINED BEHAVIOUR: Number of super-droplets in the domain cannot become larger than the "
        "size of the super-droplets' view");
    throw std::invalid_argument(err);
  }

  for (size_t kk(start); kk < start + newnsupers; ++kk) {
    totsupers(kk) = create_superdrop(gbxmaps, gbxindex);
  }
}

/*
Call to apply boundary conditions to remove and then add superdroplets to the top of the domain
abouve coord3lim.

_Note:_ totsupers is view of all superdrops (both in and out of bounds of domain).
*/
void AddSupersAtDomainTop::operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                                      const viewd_supers totsupers) const {
  const size_t ngbxs(d_gbxs.extent(0));

  bool is_supers_added = false;
  remove_superdrops_from_gridboxes(gbx);

  for (size_t ii(0); ii < ngbxs; ++ii) {  // TODO(CB) parallelise?
    auto &gbx(d_gbxs(ii));
    const auto bounds = gbxmaps.coord3bounds(gbx.get_gbxindex());
    if (bounds.second > coord3lim) {
      add_superdrops_for_gridbox(gbxmaps, gbx, totsupers);
      is_supers_added = true;
    }
  }

  if (is_supers_added) {  // resort totsupers view and set gbx references
    move_supers_between_gridboxes_again(d_gbxs, totsupers);
  }
}

/* call to create a new superdroplet for gridbox with given gbxindex */
Superdrop CreateSuperdrop::operator()(const CartesianMaps &gbxmaps, const auto gbxindex) const {
  const auto sdgbxindex = gbxindex;
  const auto coords312 = create_superdrop_coords(gbxmaps, gbxindex);
  const auto attrs = create_superdrop_attrs(gbxmaps.get_gbxvolume(gbxindex));
  const auto sdId = sdIdGen->next();

  return Superdrop(sdgbxindex, coords312[0], coords312[1], coords312[2], attrs, sdId);
}

/* create spatial coordinates for super-droplet by setting coord1 = coord2 = 0.0 and coord3 to a
random value within the gridbox's bounds */
std::array<double, 3> CreateSuperdrop::create_superdrop_coords(const CartesianMaps &gbxmaps,
                                                               const auto gbxindex) const {
  const auto bounds = gbxmaps.coord3bounds(gbxindex);
  const auto coord3 = randgen->drand(bounds.first, bounds.second);
  const auto coord1 = double{0.0 / dlc::COORD0};
  const auto coord2 = double{0.0 / dlc::COORD0};

  return std::array<double, 3>{coord3, coord1, coord2};
}

/* create attributes for a new super-droplet */
SuperdropAttrs CreateSuperdrop::create_superdrop_attrs(const double gbxvolume) const {
  const auto [xi, radius] = new_xi_radius(gbxvolume);
  const auto msol = new_msol(radius);
  const auto solute = SoluteProperties{};

  return SuperdropAttrs(solute, xi, radius, msol);
}

/* returns radius and xi for a new super-droplet by randomly sampling a distribution. */
std::pair<size_t, double> CreateSuperdrop::new_xi_radius(const double gbxvolume) const {
  const auto bin = uint64_t{randgen->urand(0, nbins)};  // index of randomly selected log10(r) bin

  const auto log10rlow = log10redges.at(bin);     // lower bound of log10(r)
  const auto log10rup = log10redges.at(bin + 1);  // upper bound of log10(r)
  const auto log10rwidth = (log10rup - log10rlow);
  const auto frac = randgen->drand(0.0, 1.0);
  const auto log10r = double{log10rlow + frac * log10rwidth};
  const auto radius = double{std::pow(10.0, log10r)};

  const auto nconc = dist.droplet_numconc_distribution(log10r, log10rup, log10rlow);
  const auto xi = (uint64_t)std::round(nconc * gbxvolume);  // cast double to uint64_t

  return std::make_pair(xi, radius);  // xi_radius
}

/* returns solute mass for a new super-droplet with a dryradius = 1nano-meter. */
double CreateSuperdrop::new_msol(const double radius) const {
  constexpr double msolconst = 4.0 * std::numbers::pi * dlc::Rho_sol / 3.0;

  if (radius < dryradius) {
    throw std::invalid_argument("new radius cannot be < dry radius of droplet");
  }

  return msolconst * dryradius * dryradius * dryradius;
}

/* returns the droplet number concentration for a bin of width log10rlow -> log10rup
from the sum of two Lognormal distributions centered on the radius at log10r. */
double TwoLognormalsDistribution::droplet_numconc_distribution(double log10r, double log10rup,
                                                               double log10rlow) const {
  const auto nconc_a = dist_a.droplet_numconc_distribution(log10r, log10rup, log10rlow);
  const auto nconc_b = dist_b.droplet_numconc_distribution(log10r, log10rup, log10rlow);
  return nconc_a + nconc_b;
}

/* normalised lognormal distribution returns the probability density of a given radius */
double LognormalDistribution::lognormal_pdf(double radius) const {
  double inverse_norm = radius * lnsigma * std::sqrt(2 * std::numbers::pi);
  double expo = std::log(radius / geomean) / lnsigma;
  return std::exp(-0.5 * expo * expo) / inverse_norm;
}

/* returns the droplet number concentration for a bin of width log10rlow -> log10rup
from a Lognormal distribution centered on the radius at log10r. */
double LognormalDistribution::droplet_numconc_distribution(double log10r, double log10rup,
                                                           double log10rlow) const {
  double delta_radius = std::pow(10.0, log10rup) - std::pow(10.0, log10rlow);
  double dnumconc_dradius = numconc * lognormal_pdf(std::pow(10.0, log10r));
  return dnumconc_dradius * delta_radius;  // number of droplets per unit volume for bin
}
