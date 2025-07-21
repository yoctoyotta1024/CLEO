/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: add_supers_at_domain_top.hpp
 * Project: movement
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * A definition of the a type satisyfing the BoundaryConditions concept
 * to use for a Cartesian Domain in MoveSupersInDomain.
 */

#ifndef LIBS_CARTESIANDOMAIN_MOVEMENT_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
#define LIBS_CARTESIANDOMAIN_MOVEMENT_ADD_SUPERS_AT_DOMAIN_TOP_HPP_

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../../cleoconstants.hpp"
#include "../../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/domainboundaries.hpp"
#include "configuration/optional_config_params.hpp"
#include "gridboxes/supersindomain.hpp"
#include "superdrops/superdrop.hpp"

namespace KCS = KokkosCleoSettings;

struct LognormalDistribution {
  double numconc; /**< number concentration of new droplets */
  double geomean; /**< geometric mean of lognormal distribution */
  double lnsigma; /**< ln(geometric sigma) of lognormal distribution */

  /* normalised lognormal distribution returns the probability density of a given radius */
  double lognormal_pdf(double radius) const;

  /* returns the droplet number concentration for a bin of width log10rlow -> log10rup
  from a Lognormal distribution centered on the radius at log10r. */
  double droplet_numconc_distribution(double log10r, double log10rup, double log10rlow) const;
};

class TwoLognormalsDistribution {
 private:
  LognormalDistribution dist_a; /**< 1st lognormal distribution for creating superdroplet xi */
  LognormalDistribution dist_b; /**< 2nd lognormal distribution for creating superdroplet xi */

 public:
  explicit TwoLognormalsDistribution(const OptionalConfigParams::AddSupersAtDomainTopParams &config)
      : dist_a({config.NUMCONC_a * dlc::VOL0, config.GEOMEAN_a / dlc::R0,
                std::log(config.geosigma_a)}),
        dist_b({config.NUMCONC_b * dlc::VOL0, config.GEOMEAN_b / dlc::R0,
                std::log(config.geosigma_b)}) {}

  /* returns the droplet number concentration for a bin of width log10rlow -> log10rup
  from the sum of two Lognormal distributions centered on the radius at log10r. */
  double droplet_numconc_distribution(double log10r, double log10rup, double log10rlow) const;
};

struct CreateSuperdrop {
 private:
  std::shared_ptr<std::mt19937> randgen; /**< pointer to random number generator */
  std::shared_ptr<Superdrop::IDType::Gen>
      sdIdGen;  /**< Pointer Superdrop::IDType object for super-droplet ID generation. */
  size_t nbins; /**< number of bins for sampling superdroplet radius */
  std::vector<double> log10redges; /**< edges of bins for superdroplet log_10(radius) */
  double dryradius;                /**< dry radius of new superdrop */
  TwoLognormalsDistribution dist;  /**< distribution for creating superdroplet xi */

  /* create spatial coordinates for super-droplet by setting coord1 = coord2 = 0.0 and coord3 to
  a random value within the gridbox's bounds */
  std::array<double, 3> create_superdrop_coords(const CartesianMaps &gbxmaps,
                                                const unsigned int gbxindex) const;

  /* create attributes for a new super-droplet */
  SuperdropAttrs create_superdrop_attrs(const double gbxvolume) const;

  /* returns radius and xi for a new super-droplet by randomly sampling a distribution. */
  std::pair<size_t, double> new_xi_radius(const double gbxvolume) const;

  /* returns solute mass for a new super-droplet with a dryradius = 1nano-meter. */
  double new_msol(const double radius) const;

 public:
  /* call to create a new superdroplet for gridbox with given gbxindex */
  explicit CreateSuperdrop(const OptionalConfigParams::AddSupersAtDomainTopParams &config);

  Superdrop operator()(const CartesianMaps &gbxmaps, const unsigned int gbxindex) const;
};

/*
 * struct satisfying BoundaryConditions concept for applying domain boundary conditions which
 * add superdroplets to gridboxes above a certain height.
 */
struct AddSupersAtDomainTop {
 private:
  size_t newnsupers; /**< number of superdroplets to add to gridboxes above coord3lim */
  double coord3lim;  /**< gridboxes with upper bound > coord3lim get new super-droplets */
  CreateSuperdrop create_superdrop; /**< methods to create a new superdrop */

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
  */
  SupersInDomain apply(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                       SupersInDomain &allsupers) const;
};

#endif  // LIBS_CARTESIANDOMAIN_MOVEMENT_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
