/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: coalescence.hpp
 * Project: collisions
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Class and function to enact collision-coalescence events in the super-droplet model
 * according to Shima et al. 2009. DoCoalescence struct satisfies PairEnactX concept
 * used by DoCollisions.
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_COALESCENCE_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_COALESCENCE_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <cstdint>
#include <functional>

#include "../microphysicalprocess.hpp"
#include "../superdrop.hpp"
#include "./collisions.hpp"

/**
 * @brief Raises an error if the multiplicity of the super-droplet is 0.
 *
 * This function checks if the given superdrop is null by verifying its multiplicity.
 * It asserts that the multiplicity (`xi`) of the superdrop is greater than 0.
 * If the multiplicity is 0, an assertion error is raised.
 *
 * @param drop A reference to the `Superdrop` object to be checked.
 * @return Always returns 0 (=false).
 */
KOKKOS_INLINE_FUNCTION bool is_null_superdrop(const Superdrop &drop) {
  assert((drop.get_xi() > 0) && "superdrop xi < 1, null drop in coalescence");
  return 0;
}

struct DoCoalescence {
 private:
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
  KOKKOS_FUNCTION void twin_superdroplet_coalescence(const uint64_t gamma, Superdrop &drop1,
                                                     Superdrop &drop2) const;

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
  KOKKOS_FUNCTION void different_superdroplet_coalescence(const uint64_t gamma, Superdrop &drop1,
                                                          Superdrop &drop2) const;

 public:
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
  KOKKOS_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2, const double prob, const double phi) const;

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
  KOKKOS_FUNCTION uint64_t coalescence_gamma(const uint64_t xi1, const uint64_t xi2,
                                             const double prob, const double phi) const;

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
  KOKKOS_FUNCTION bool coalesce_superdroplet_pair(const uint64_t gamma, Superdrop &drop1,
                                                  Superdrop &drop2) const;
};

/**
 * @brief Constructs a microphysical process for collision-coalescence of superdroplets.
 *
 * This function constructs a microphysical process for collision-coalescence of superdroplets with
 * a constant timestep and probability of collision-coalescence determined by 'collcoalprob'.
 *
 * @tparam Probability Type satisfying the PairProbability concept.
 * @param interval The constant timestep interval.
 * @param int2realtime A function that converts an integer timestep to real time.
 * @param collcoalprob The probability of collision-coalescence.
 * @return An instance of MicrophysicalProcess for collision-coalescence.
 */
template <PairProbability Probability>
inline MicrophysicalProcess auto CollCoal(const unsigned int interval,
                                          const std::function<double(unsigned int)> int2realtime,
                                          const Probability collcoalprob) {
  const auto DELT = int2realtime(interval);

  const DoCoalescence coal{};
  const MicrophysicsFunc auto colls =
      DoCollisions<Probability, DoCoalescence>(DELT, collcoalprob, coal);

  return ConstTstepMicrophysics(interval, colls);
}

#endif  // LIBS_SUPERDROPS_COLLISIONS_COALESCENCE_HPP_
