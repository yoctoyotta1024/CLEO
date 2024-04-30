/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: withreset.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 19th December 2023
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
 * Motion of a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * the sdgbxindex updated accordingly for a
 * cartesian domain with finite/periodic boundary
 * conditions and reset of superdroplets that leave
 * the domain through coord3 domain boundaries
 */

// TODO(CB) Delete file (!)

#ifndef LIBS_CARTESIANDOMAIN_WITHRESET_HPP_
#define LIBS_CARTESIANDOMAIN_WITHRESET_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <functional>
#include <random>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/cartesianmotion.hpp"
#include "gridboxes/predcorrmotion.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"
#include "superdrops/urbg.hpp"

namespace dlc = dimless_constants;

struct ProbDistrib {
 private:
  double REFF;
  double nueff;
  double n0const;

  /* returns normalised probability density distribution,
  ie. probability of radius in range r -> r+ dr, such that
  integral over all radii = 1. Distribution is gamma
  distribution for cloud droplets using parameters
  from Poertge et al. 2023 for shallow cumuli (figure 12), ie.
  with typical values: reff = 7e-6 m, and nueff = 0.08.
  RADIUS has dimensions [m] */
  KOKKOS_FUNCTION double probdens_distrib(const double RADIUS) const {
    const auto term1 = double{Kokkos::pow(RADIUS, ((1.0 - 3.0 * nueff) / nueff))};
    const auto term2 = double{Kokkos::exp(-RADIUS / (REFF * nueff))};

    const auto probdens = double{n0const * term1 * term2};  // dn_dr [prob m^-1]

    return probdens;  // normalised probability in range r -> r+dr
  }

 public:
  ProbDistrib() : REFF(7e-6), nueff(0.08), n0const(0.0) {
    const auto xp = double{(1.0 - 2.0 * nueff) / nueff};
    const auto valxp = double{Kokkos::pow(REFF * nueff, -xp)};
    n0const = valxp / Kokkos::tgamma(xp);
  }

  /* returns probability of radius in range r -> r+ dr, such that
  integral of probability density dsitribution over all radii = 1 */
  KOKKOS_FUNCTION double operator()(const double radius, const double rlow,
                                    const double rup) const {
    const auto RADIUS = radius * dlc::R0;        // dimensionalised radius [m]
    const auto DELTAR = (rup - rlow) * dlc::R0;  // dimensionalised bin wdith [m]

    const auto prob = probdens_distrib(RADIUS) * DELTAR;

    return prob;  // probability of radius
  }
};

struct ResetSuperdrop {
  GenRandomPool genpool4reset;
  Kokkos::View<double[101]> log10redges;  // edges to radius bins
  Kokkos::pair<size_t, size_t> gbxidxs;
  uint64_t nbins;
  ProbDistrib prob_distrib;

  ResetSuperdrop(const size_t ngbxs, const size_t ngbxs4reset)
      : genpool4reset(std::random_device {}()),
        log10redges("log10redges"),
        gbxidxs({ngbxs - ngbxs4reset, ngbxs}),
        nbins(log10redges.extent(0) - 1),
        prob_distrib(ProbDistrib()) {
    /* make redges linearly spaced in log10(R) space */
    auto h_log10redges = Kokkos::create_mirror_view(log10redges);
    const auto log10rmin = double{Kokkos::log10(5e-6 / dlc::R0)};    // lowest edge of radius bins
    const auto log10rmax = double{Kokkos::log10(1.5e-5 / dlc::R0)};  // highest edge of radius bins
    const auto log10deltar = double{(log10rmax - log10rmin) / nbins};
    for (size_t i(0); i < nbins + 1; ++i) {
      h_log10redges(i) = log10rmin + i * log10deltar;
    }
    Kokkos::deep_copy(log10redges, h_log10redges);
  }

  /* randomly update position of superdroplet by
  randomly selecting a gbxindex from gbxidxs and then
  randomly selecting a coord3 with that gbx's bounds */
  KOKKOS_FUNCTION unsigned int reset_position(const CartesianMaps &gbxmaps, URBG<ExecSpace> &urbg,
                                              Superdrop &drop) const {
    const auto sdgbxindex =
        urbg(gbxidxs.first,
             gbxidxs.second);  // randomly selected gbxindex in range {incl., excl.}

    const auto bounds = gbxmaps.coord3bounds(sdgbxindex);
    const auto coord3 = urbg.drand(bounds.first, bounds.second);  // random coord within gbx bounds

    drop.set_sdgbxindex(sdgbxindex);
    drop.set_coord3(coord3);

    return sdgbxindex;
  }

  /* reset radius and multiplicity of superdroplet
  by randomly sampling from binned distributions */
  KOKKOS_FUNCTION void reset_attributes(const double gbxvol, URBG<ExecSpace> &urbg,
                                        Superdrop &drop) const {
    const auto bin = uint64_t{urbg(0, nbins)};   // index of randomly selected log10(r) bin
    const auto log10rlow = log10redges(bin);     // lower bound of log10(r)
    const auto log10rup = log10redges(bin + 1);  // upper bound of log10(r)

    const auto radius = new_radius(log10rlow, log10rup, urbg);
    const auto xi = new_xi(gbxvol, log10rlow, log10rup, radius);
    const auto msol = new_msol(radius);

    drop.set_msol(msol);
    drop.set_radius(radius * 1.00000001);  // radius 1e-6 % larger than sampled dryradius
    drop.set_xi(xi);
  }

  /* returns msol given dry radius */
  KOKKOS_FUNCTION double new_msol(const double dryradius) const {
    constexpr double msolconst = 4.0 * Kokkos::numbers::pi * dlc::Rho_sol / 3.0;

    return msolconst * dryradius * dryradius * dryradius;
  }

  /* returns radius from within bin of uniform
  distiribution in log10(r) space */
  KOKKOS_FUNCTION double new_radius(const double log10rlow, const double log10rup,
                                    URBG<ExecSpace> &urbg) const {
    const auto frac = urbg.drand(0.0, 1.0);
    const auto log10r = double{log10rlow + frac * (log10rup - log10rlow)};
    const auto radius = double{Kokkos::pow(10.0, log10r)};

    return radius;
  }

  /* returns xi given value of normalised probability
  distribution at radius and the bin width */
  KOKKOS_FUNCTION uint64_t new_xi(const double gbxvol, const double log10rlow,
                                  const double log10rup, const double radius) const {
    constexpr double numconc = 100000000 * dlc::VOL0;  // 100/cm^3, non-dimensionalised

    const auto rlow = double{Kokkos::pow(10.0, log10rlow)};
    const auto rup = double{Kokkos::pow(10.0, log10rup)};

    const auto prob = prob_distrib(radius, rlow, rup);
    const auto xi = double{prob * numconc * gbxvol};

    return (uint64_t)Kokkos::round(xi);
  }

  KOKKOS_FUNCTION unsigned int operator()(const CartesianMaps &gbxmaps, Superdrop &drop) const {
    URBG<ExecSpace> urbg{genpool4reset.get_state()};  // thread safe random number generator

    const auto sdgbxindex = reset_position(gbxmaps, urbg, drop);
    const auto gbxvol = gbxmaps.get_gbxvolume(sdgbxindex);
    reset_attributes(gbxvol, urbg, drop);

    genpool4reset.free_state(urbg.gen);

    return sdgbxindex;
  }
};

KOKKOS_FUNCTION unsigned int change_if_coord3nghbr_withreset(const ResetSuperdrop &reset_superdrop,
                                                             const CartesianMaps &gbxmaps,
                                                             unsigned int idx, Superdrop &drop);

/* wrapper of functions for use in PredCorrMotion's
ChangeToNghbr type for deciding if a superdroplet should move
to a neighbouring gbx in a cartesian domain and then updating the
superdroplet appropriately. Struct has three functions, one
for each direction (coord3 = z, coord1 = x, coord2 = y). For each,
the superdrop's coord is compared to gridbox bounds given by gbxmaps
for the current gbxindex 'idx'. If superdrop coord lies outside
bounds, forward or backward neighbour functions are called to
update sdgbxindex (and possibly other superdrop attributes).
Struct is same as CartesianChangeIfNghbr except for in
coord3(...){...} function */
struct CartesianChangeIfNghbrWithReset {
  ResetSuperdrop reset_superdrop;

  CartesianChangeIfNghbrWithReset(const size_t ngbxs, const size_t ngbxs4reset)
      : reset_superdrop(ResetSuperdrop(ngbxs, ngbxs4reset)) {}

  KOKKOS_INLINE_FUNCTION unsigned int coord3(const CartesianMaps &gbxmaps, unsigned int idx,
                                             Superdrop &drop) const {
    return change_if_coord3nghbr_withreset(reset_superdrop, gbxmaps, idx, drop);
  }

  KOKKOS_INLINE_FUNCTION unsigned int coord1(const CartesianMaps &gbxmaps, unsigned int idx,
                                             Superdrop &drop) const {
    return change_if_coord1nghbr(gbxmaps, idx, drop);
  }

  KOKKOS_INLINE_FUNCTION unsigned int coord2(const CartesianMaps &gbxmaps, unsigned int idx,
                                             Superdrop &drop) const {
    return change_if_coord2nghbr(gbxmaps, idx, drop);
  }
};

/* returned type satisfies motion concept for motion of a
superdroplet using a predictor-corrector method to update
a superdroplet's coordinates and then updating it's
sdgbxindex as appropriate for a cartesian domain */
template <VelocityFormula TV>
inline PredCorrMotion<CartesianMaps, TV, CartesianChangeIfNghbrWithReset, CartesianCheckBounds>
CartesianMotionWithReset(const unsigned int motionstep,
                         const std::function<double(unsigned int)> int2time, const TV terminalv,
                         const size_t ngbxs, const size_t ngbxs4reset) {
  const auto cin = CartesianChangeIfNghbrWithReset(ngbxs, ngbxs4reset);
  return PredCorrMotion<CartesianMaps, TV, CartesianChangeIfNghbrWithReset, CartesianCheckBounds>(
      motionstep, int2time, terminalv, cin, CartesianCheckBounds{});
}

#endif  // LIBS_CARTESIANDOMAIN_WITHRESET_HPP_
