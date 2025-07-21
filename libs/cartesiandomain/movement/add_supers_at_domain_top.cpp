/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: add_supers_at_domain_top.cpp
 * Project: movement
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Implementation of AddSupersAtDomainTop type satisyfing the BoundaryConditions concept
 * to use for a Cartesian Domain in MoveSupersInDomain.
 */

#include "./add_supers_at_domain_top.hpp"

Kokkos::View<unsigned int *> remove_superdrops_from_gridboxes(const CartesianMaps &gbxmaps,
                                                              const viewd_gbx d_gbxs,
                                                              const subviewd_supers domainsupers,
                                                              const double coord3lim);
viewd_supers create_newsupers_for_gridboxes(const CartesianMaps &gbxmaps,
                                            const CreateSuperdrop &create_superdrop,
                                            Kokkos::View<unsigned int *> gbxindexes,
                                            const size_t newnsupers_pergbx);
void add_superdrops_for_gridboxes(const SupersInDomain &allsupers,
                                  const viewd_constsupers newsupers);
SupersInDomain move_supers_between_gridboxes_again(const viewd_gbx d_gbxs,
                                                   SupersInDomain &allsupers);

/*
Call to apply boundary conditions to remove and then add superdroplets to the top of the domain
above coord3lim.
*/
SupersInDomain AddSupersAtDomainTop::apply(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                                           SupersInDomain &allsupers) const {
  const auto gbxindexes_for_newsupers =
      remove_superdrops_from_gridboxes(gbxmaps, d_gbxs, allsupers.domain_supers(), coord3lim);

  auto newsupers_for_gridboxes = create_newsupers_for_gridboxes(
      gbxmaps, create_superdrop, gbxindexes_for_newsupers, newnsupers);

  add_superdrops_for_gridboxes(allsupers, newsupers_for_gridboxes);

  allsupers = move_supers_between_gridboxes_again(d_gbxs, allsupers);

  return allsupers;
}

/* (re)sorting supers based on their gbxindexes and then updating the span (gbx refs) for each
gridbox accordingly.

_Note:_ totsupers is view of all superdrops (both in and out of bounds of domain).

Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}.
*/
SupersInDomain move_supers_between_gridboxes_again(const viewd_gbx d_gbxs,
                                                   SupersInDomain &allsupers) {
  allsupers.sort_totsupers(d_gbxs);

  const auto domainsupers = allsupers.domain_supers();
  const size_t ngbxs(d_gbxs.extent(0));
  Kokkos::parallel_for(
      "move_supers_between_gridboxes_again", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
      KOKKOS_LAMBDA(const size_t ii) { d_gbxs(ii).supersingbx.set_refs(domainsupers); });

  return allsupers;
}

/* set super-droplet sdgbxindex to out of bounds value if superdrop coord3 > coord3lim.
Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t kk(0); kk < supers.extent(0); ++kk){[...]}.
*/
KOKKOS_FUNCTION
void remove_superdrop_above_coord3lim(const TeamMember &team_member,
                                      const subviewd_supers domainsupers, const Gridbox &gbx,
                                      const double coord3lim) {
  const auto supers = gbx.supersingbx(domainsupers);
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team_member, supers.extent(0)), [supers, coord3lim](const size_t kk) {
        if (supers(kk).get_coord3() >= coord3lim) {
          supers(kk).set_sdgbxindex(LIMITVALUES::oob_gbxindex);  // remove super-droplet from domain
        }
      });
}

/* for gridboxes with coordinates above coord3lim, set super-droplet sdgbxindex to out of bounds
value if superdrop coord3 > coord3lim.

Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}.

Function returns view of all the gridboxes indexes in d_gbxs where the value of the gridbox index
has been replaced by outofbounds_gbxindex unless superdrops were removed from that gridbox.
*/
Kokkos::View<unsigned int *> remove_superdrops_from_gridboxes(const CartesianMaps &gbxmaps,
                                                              const viewd_gbx d_gbxs,
                                                              const subviewd_supers domainsupers,
                                                              const double coord3lim) {
  const size_t ngbxs(d_gbxs.extent(0));
  auto gbxindexes_of_removedsupers = Kokkos::View<unsigned int *>("gbxs_of_removedsupers", ngbxs);
  Kokkos::parallel_for(
      "remove_superdrops", TeamPolicy(ngbxs, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const TeamMember &team_member) {
        const auto ii = team_member.league_rank();

        const auto ubound = gbxmaps.coord3bounds(d_gbxs(ii).get_gbxindex()).second;
        if (ubound > coord3lim) {
          remove_superdrop_above_coord3lim(team_member, domainsupers, d_gbxs(ii), coord3lim);
          gbxindexes_of_removedsupers(ii) = d_gbxs(ii).get_gbxindex();  // add newsupers
        } else {
          gbxindexes_of_removedsupers(ii) = LIMITVALUES::oob_gbxindex;  // don't add newsupers
        }
      });

  return gbxindexes_of_removedsupers;
}

/* Given a view of gridboxes where the value of the gridbox index has been replaced by
outofbounds_gbxindex unless superdrops should be added to that gridbox,
count the total number of new superdroplets to create.

Kokkos::parallel_reduce([...]) is equivalent in serial to suming up newnsupers_total in loop:
for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}.
*/
size_t total_newnsupers_to_create(Kokkos::View<unsigned int *> gbxindexes,
                                  const size_t newnsupers_pergbx) {
  auto newnsupers_total = size_t{0};
  Kokkos::parallel_reduce(
      "newnsupers_total", Kokkos::RangePolicy<ExecSpace>(0, gbxindexes.extent(0)),
      KOKKOS_LAMBDA(const size_t ii, size_t &nsupers) {
        if (gbxindexes(ii) != LIMITVALUES::oob_gbxindex) {
          nsupers = newnsupers_pergbx;
        } else {
          nsupers = 0;
        }
      },
      newnsupers_total);

  return newnsupers_total;
}

/* Given a view of gridboxes where the value of the gridbox index has been replaced by
outofbounds_gbxindex unless superdrops should be added to that gridbox, function create 'newnsupers'
new superdroplets per gridbox by calling the create_superdrop function on host and then
copies resultant view to device memory.
*/
viewd_supers create_newsupers_for_gridboxes(const CartesianMaps &gbxmaps,
                                            const CreateSuperdrop &create_superdrop,
                                            Kokkos::View<unsigned int *> gbxindexes,
                                            const size_t newnsupers_pergbx) {
  auto newsupers =
      viewd_supers("newsupers", total_newnsupers_to_create(gbxindexes, newnsupers_pergbx));
  auto h_newsupers = Kokkos::create_mirror_view(newsupers);

  auto h_gbxindexes = Kokkos::create_mirror_view(gbxindexes);
  Kokkos::deep_copy(h_gbxindexes, gbxindexes);

  auto nn = size_t{0};  // number of super_droplets created
  for (size_t ii(0); ii < h_gbxindexes.extent(0); ++ii) {
    if (h_gbxindexes(ii) != LIMITVALUES::oob_gbxindex) {
      for (size_t kk(0); kk < newnsupers_pergbx; ++kk) {
        h_newsupers(nn) = create_superdrop(gbxmaps, h_gbxindexes(ii));
        ++nn;
      }
    }
  }
  Kokkos::deep_copy(newsupers, h_newsupers);

  assert((newsupers.extent(0) == nn) &&
         "total number of superdrops created must equal newsupers view size");

  return newsupers;
}

/* returns copy of 1 gridbox within a view on host memory of the ii'th gridbox in a
device view 'd_gbxs' */
viewh_constgbx hostcopy_one_gridbox(const viewd_constgbx d_gbxs, const size_t ii) {
  const auto d_gbx = viewd_gbx("gbx", 1);
  Kokkos::parallel_for(
      "copy_gbx", 1, KOKKOS_LAMBDA(const unsigned int i) { d_gbx(0) = d_gbxs(ii); });

  auto h_gbx = Kokkos::create_mirror_view(d_gbx);
  Kokkos::deep_copy(h_gbx, d_gbx);
  return h_gbx;
}

/* throws error if the size of newnsupers + oldnsupers > total space in totsupers view */
size_t check_space_in_totsupers(const SupersInDomain &allsupers,
                                const viewd_constsupers newsupers) {
  const auto totsupers = allsupers.get_totsupers_readonly();
  const auto oldnsupers = allsupers.domain_nsupers();
  if (oldnsupers + newsupers.extent(0) > totsupers.extent(0)) {
    const auto err = std::string(
        "UNDEFINED BEHAVIOUR: Number of super-droplets in the domain cannot become larger than the "
        "size of the super-droplets' view");
    throw std::invalid_argument(err);
  }

  return oldnsupers;
}

/* check there is space in totsupers for newsupers, then append superdrops in newsupers to end of
totsupers view */
void add_superdrops_for_gridboxes(const SupersInDomain &allsupers,
                                  const viewd_constsupers newsupers) {
  const auto totsupers = allsupers.get_totsupers();
  const auto og_ntotsupers = check_space_in_totsupers(allsupers, newsupers);

  Kokkos::parallel_for(
      "add_superdrops", Kokkos::RangePolicy<ExecSpace>(0, newsupers.extent(0)),
      KOKKOS_LAMBDA(const size_t kk) { totsupers(kk + og_ntotsupers) = newsupers(kk); });
}

/* returns host copy of {lower, upper} coord3 boundaries from gbxmaps for 'gbxindex' on device */
Kokkos::pair<double, double> hostcopy_coord3bounds(const CartesianMaps &gbxmaps,
                                                   const unsigned int gbxindex) {
  const Kokkos::View<Kokkos::pair<double, double>[1]> d_bound("d_bound");
  Kokkos::parallel_for(
      "copy_coord3bound", 1,
      KOKKOS_LAMBDA(const unsigned int i) { d_bound(0) = gbxmaps.coord3bounds(gbxindex); });

  auto h_bound = Kokkos::create_mirror_view(d_bound);
  Kokkos::deep_copy(h_bound, d_bound);
  return h_bound(0);
}

/* call to create a new superdroplet for gridbox with given gbxindex */
CreateSuperdrop::CreateSuperdrop(const OptionalConfigParams::AddSupersAtDomainTopParams &config)
    : randgen(std::make_shared<std::mt19937>(std::random_device {}())),
      sdIdGen(std::make_shared<Superdrop::IDType::Gen>(config.initnsupers)),
      nbins(config.newnsupers),
      log10redges(),
      dryradius(config.DRYRADIUS / dlc::R0),
      dist(config) {
  const auto log10rmin = std::log10(config.MINRADIUS / dlc::R0);
  const auto log10rmax = std::log10(config.MAXRADIUS / dlc::R0);
  const auto log10deltar = double{(log10rmax - log10rmin) / nbins};
  for (size_t nn(0); nn < nbins + 1; ++nn) {
    log10redges.push_back(log10rmin + nn * log10deltar);
  }
}

/* call to create a new superdroplet for gridbox with given gbxindex */
Superdrop CreateSuperdrop::operator()(const CartesianMaps &gbxmaps,
                                      const unsigned int gbxindex) const {
  const auto sdgbxindex = gbxindex;
  const auto coords312 = create_superdrop_coords(gbxmaps, gbxindex);
  const auto attrs = create_superdrop_attrs(gbxmaps.get_gbxvolume(gbxindex));
  const auto sdId = sdIdGen->next();

  return Superdrop(sdgbxindex, coords312[0], coords312[1], coords312[2], attrs, sdId);
}

/* create spatial coordinates for super-droplet by setting coord1 = coord2 = 0.0 and coord3 to a
random value within the gridbox's bounds */
std::array<double, 3> CreateSuperdrop::create_superdrop_coords(const CartesianMaps &gbxmaps,
                                                               const unsigned int gbxindex) const {
  const auto bounds = hostcopy_coord3bounds(gbxmaps, gbxindex);
  auto dist = std::uniform_real_distribution<double>(bounds.first, bounds.second);
  const double coord3 = dist(*randgen);

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
  auto uintdist = std::uniform_int_distribution<uint64_t>(0, nbins - 1);
  const uint64_t bin = uintdist(*randgen);  // index of randomly selected log10(r) bin

  const auto log10rlow = log10redges.at(bin);     // lower bound of log10(r)
  const auto log10rup = log10redges.at(bin + 1);  // upper bound of log10(r)
  const auto log10rwidth = (log10rup - log10rlow);
  auto dbldist = std::uniform_real_distribution<double>(0.0, 1.0);
  const auto frac = dbldist(*randgen);

  const auto log10r = double{log10rlow + frac * log10rwidth};
  const auto radius = double{std::pow(10.0, log10r)};

  const auto nconc = dist.droplet_numconc_distribution(log10r, log10rup, log10rlow);
  const auto xi = static_cast<uint64_t>(std::max(std::round(nconc * gbxvolume), 1.0));

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
