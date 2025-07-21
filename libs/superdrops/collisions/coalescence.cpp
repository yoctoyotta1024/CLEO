/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: coalescence.cpp
 * Project: collisions
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functionality to enact collision-coalescence events in the super-droplet model
 * according to Shima et al. 2009. DoCoalescence struct satisfies PairEnactX concept
 * used by DoCollisions.
 */

#include "./coalescence.hpp"

/**
 * @brief Operator used as an adaptor such that DoCoalescence satisfies the PairEnactX concept
 * and so can be used as the EnactCollision function-like object in the DoCollisions struct.
 *
 * This operator calls functions to enact the collision-coalescence of two super-droplets.
 *
 * @param drop1 The first super-droplet.
 * @param drop2 The second super-droplet.
 * @param prob The probability of collision-coalescence.
 * @param phi Random number in the range [0.0, 1.0].
 * @return boolean=true if collision-coalescence resulted in null superdrops.
 */
KOKKOS_FUNCTION bool DoCoalescence::operator()(Superdrop &drop1, Superdrop &drop2,
                                               const double prob, const double phi) const {
  /* 1. calculate gamma factor for collision-coalescence  */
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();
  const auto gamma = coalescence_gamma(xi1, xi2, prob, phi);

  /* 2. enact collision-coalescence on pair
  of superdroplets if gamma is not zero */
  if (gamma != 0) {
    return coalesce_superdroplet_pair(gamma, drop1, drop2);
  }

  return 0;
}

/**
 * @brief Calculates the value of the gamma factor in Monte Carlo collision-coalescence.
 *
 * This function calculates the value of the gamma factor used in Monte Carlo
 * collision-coalescence as described in Shima et al. 2009.
 *
 * @param xi1 The multiplicity of the first super-droplet.
 * @param xi2 The multiplicity of the second super-droplet.
 * @param prob The probability of collision-coalescence.
 * @param phi Random number in the range [0.0, 1.0].
 * @return The calculated value of the coalescence gamma factor.
 */
KOKKOS_FUNCTION uint64_t DoCoalescence::coalescence_gamma(const uint64_t xi1, const uint64_t xi2,
                                                          const double prob,
                                                          const double phi) const {
  uint64_t gamma = Kokkos::floor(prob);  // if phi >= (prob - floor(prob))
  if (phi < (prob - gamma)) {
    ++gamma;
  }

  const auto maxgamma = xi1 / xi2;  // same as floor() for positive ints

  return Kokkos::fmin(gamma, maxgamma);
}

/**
 * @brief Coalesces a pair of superdroplets.
 *
 * This function coalesces a pair of superdroplets by changing their multiplicity, radius, and
 * solute mass according to Shima et al. 2009 Section 5.1.3. part (5).
 *
 * @param gamma The coalescence gamma factor.
 * @param drop1 The first superdroplet.
 * @param drop2 The second superdroplet.
 * @return True if coalescence results in a null superdroplet, false otherwise.
 */
KOKKOS_FUNCTION bool DoCoalescence::coalesce_superdroplet_pair(const uint64_t gamma,
                                                               Superdrop &drop1,
                                                               Superdrop &drop2) const {
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();

  if (xi1 - gamma * xi2 > 0) {
    different_superdroplet_coalescence(gamma, drop1, drop2);
    return 0;
  } else if (xi1 - gamma * xi2 == 0) {
    twin_superdroplet_coalescence(gamma, drop1, drop2);

    /* if xi1 = xi2 = 1 before coalesence, then xi1=0 now */
    return is_null_superdrop(drop1);
  }

  assert((xi1 >= gamma * xi2) &&
         "something undefined occured "
         "during colllision-coalescence");
  return 0;
}

/**
 * @brief Enacts coalescence of a pair of superdroplets where xi1 = gamma*xi2.
 *
 * This function coalesces a pair of superdroplets where drop1.get_xi() = gamma* drop2.get_xi() by
 * making twin superdroplets with the same xi, radius, and solute mass. It implements Shima et al.
 * 2009 Section 5.1.3. part (5) option (b).
 *
 * In rare case where xi1 = xi2 = gamma = 1, new_xi of drop1 = 0 and drop1 should be removed from
 * domain.
 *
 * _Note:_ Implicit casting of gamma (and therefore droplets' xi values) from uint64_t to double.
 *
 * @param gamma The coalescence gamma factor.
 * @param drop1 The first superdroplet.
 * @param drop2 The second superdroplet.
 */
KOKKOS_FUNCTION void DoCoalescence::twin_superdroplet_coalescence(const uint64_t gamma,
                                                                  Superdrop &drop1,
                                                                  Superdrop &drop2) const {
  assert((drop1.get_xi() == gamma * drop2.get_xi()) && "condition for twin coalescence not met");

  const auto old_xi = drop2.get_xi();  // = drop1.xi
  const auto new_xi = old_xi / 2;      // same as floor() for positive ints

  assert((new_xi < old_xi) && "coalescence must decrease multiplicity");

  const auto new_rcubed = double{drop2.rcubed() + gamma * drop1.rcubed()};
  const auto new_r = double{Kokkos::pow(new_rcubed, (1.0 / 3.0))};

  const auto new_msol = double{drop2.get_msol() + gamma * drop1.get_msol()};

  drop1.set_xi(new_xi);
  drop2.set_xi(old_xi - new_xi);

  drop1.set_radius(new_r);
  drop2.set_radius(new_r);

  drop1.set_msol(new_msol);
  drop2.set_msol(new_msol);
}

/**
 * @brief Coalesces a pair of superdroplets where xi1 > gamma*xi2.
 *
 * This function coalesces a pair of superdroplets where xi1 > gamma*xi2 by growing the radius and
 * mass of drop2 via decreasing the multiplicity of drop1. It implements Shima et al. 2009
 * Section 5.1.3. part (5) option (a).
 *
 * _Note:_ Implicit casting of gamma (i.e. therefore droplets' xi values) from uint64_t to double.
 *
 * @param gamma The coalescence gamma factor.
 * @param drop1 The first superdroplet.
 * @param drop2 The second superdroplet.
 */
KOKKOS_FUNCTION void DoCoalescence::different_superdroplet_coalescence(const uint64_t gamma,
                                                                       Superdrop &drop1,
                                                                       Superdrop &drop2) const {
  assert((drop1.get_xi() > gamma * drop2.get_xi()) && "condition on xis for coalescence not met");

  const auto new_xi = drop1.get_xi() - gamma * drop2.get_xi();

  assert((new_xi < drop1.get_xi()) && "coalescence must decrease multiplicity");

  drop1.set_xi(new_xi);

  const auto new_rcubed = double{drop2.rcubed() + gamma * drop1.rcubed()};

  drop2.set_radius(Kokkos::pow(new_rcubed, (1.0 / 3.0)));
  drop2.set_msol(drop2.get_msol() + gamma * drop1.get_msol());
}
