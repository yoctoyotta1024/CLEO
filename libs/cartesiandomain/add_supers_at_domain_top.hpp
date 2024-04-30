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
 * Last Modified: Tuesday 30th April 2024
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
#include <Kokkos_Random.hpp>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/domainboundaries.hpp"
#include "gridboxes/sortsupers.hpp"
#include "initialise/optional_config_params.hpp"
#include "superdrops/superdrop.hpp"

class LognormalDistribution {
 private:
  double numconc; /**< number concentration of new droplets */
  double geomean; /**< geometric mean of lognormal distribution */
  double lnsigma; /**< ln(geometric sigma) of lognormal distribution */

  /* normalised lognormal distribution returns the probability density of a given radius */
  double lognormal_pdf(double radius) const;

 public:
  explicit LognormalDistribution(const OptionalConfigParams::AddSupersAtDomainTopParams &config)
      : numconc(config.NUMCONC * dlc::VOL0),
        geomean(config.GEOMEAN / dlc::R0),
        lnsigma(std::log(config.geosigma)) {}

  /* returns the droplet number concentration for a bin of width log10rlow -> log10rup
  from a Lognormal distribution centered on the radius at log10r. */
  double droplet_numconc_distribution(double log10r, double log10rup, double log10rlow) const;
};

struct CreateSuperdrop {
 private:
  std::shared_ptr<Kokkos::Random_XorShift64<HostSpace>>
      randgen; /**< pointer to Kokkos random number generator */
  std::shared_ptr<Superdrop::IDType::Gen>
      sdIdGen;  /**< Pointer Superdrop::IDType object for super-droplet ID generation. */
  size_t nbins; /**< number of bins for sampling superdroplet radius */
  std::vector<double> log10redges; /**< edges of bins for superdroplet log_10(radius) */
  double dryradius;                /**< dry radius of new superdrop */
  LognormalDistribution lndist;    /**< lognormal distribution for creating superdroplet xi */

  /* create spatial coordinates for super-droplet by setting coord1 = coord2 = 0.0 and coord3 to a
  random value within the gridbox's bounds */
  std::array<double, 3> create_superdrop_coords(const CartesianMaps &gbxmaps,
                                                const auto gbxindex) const;

  /* create attributes for a new super-droplet */
  SuperdropAttrs create_superdrop_attrs(const double gbxvolume) const;

  /* returns radius and xi for a new super-droplet by randomly sampling a distribution. */
  std::pair<size_t, double> new_xi_radius(const double gbxvolume) const;

  /* returns solute mass for a new super-droplet with a dryradius = 1nano-meter. */
  double new_msol(const double radius) const;

 public:
  /* call to create a new superdroplet for gridbox with given gbxindex */
  explicit CreateSuperdrop(const OptionalConfigParams::AddSupersAtDomainTopParams &config)
      : randgen(std::make_shared<Kokkos::Random_XorShift64<HostSpace>>(std::random_device {}())),
        sdIdGen(std::make_shared<Superdrop::IDType::Gen>(config.initnsupers)),
        nbins(config.newnsupers),
        log10redges(),
        dryradius(config.DRYRADIUS / dlc::R0),
        lndist(config) {
    const auto log10rmin = std::log10(config.MINRADIUS / dlc::R0);
    const auto log10rmax = std::log10(config.MAXRADIUS / dlc::R0);
    const auto log10deltar = double{(log10rmax - log10rmin) / nbins};
    for (size_t nn(0); nn < nbins + 1; ++nn) {
      log10redges.push_back(log10rmin + nn * log10deltar);
    }
  }

  Superdrop operator()(const CartesianMaps &gbxmaps, const auto gbxindex) const;
};

struct AddSupersAtDomainTop {
 private:
  size_t newnsupers; /**< number of superdroplets to add to gridboxes above coord3lim */
  double coord3lim;  /**< gridboxes with upper bound > coord3lim get new super-droplets */
  CreateSuperdrop create_superdrop; /**< methods to create a new superdrop */

  /* set super-droplet sdgbxindex to out of bounds value if superdrop coord3 > coord3lim */
  void remove_superdrops_from_gridbox(const Gridbox &gbx) const;

  /* create 'newnsupers' number of new superdroplets from the create_superdrop function */
  void add_superdrops_for_gridbox(const CartesianMaps &gbxmaps, const Gridbox &gbx,
                                  const viewd_supers totsupers) const;

 public:
  /* New super-droplets are added to domain with coord3 >= COORD3LIM [m]. Note generation of
   * nextsdId assumes it is the only method creating super-droplets - otherwise created sdId may not
   * be unique*/
  explicit AddSupersAtDomainTop(const OptionalConfigParams::AddSupersAtDomainTopParams &config)
      : newnsupers(config.newnsupers),
        coord3lim(config.COORD3LIM / dlc::COORD0),
        create_superdrop(config) {}

  /*
  Call to apply boundary conditions to remove and then add superdroplets to the top of the domain
  abouve coord3lim.

  _Note:_ totsupers is view of all superdrops (both in and out of bounds of domain).
  */
  void operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                  const viewd_supers totsupers) const;
};

#endif  // LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
