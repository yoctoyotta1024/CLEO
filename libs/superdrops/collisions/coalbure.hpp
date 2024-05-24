/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: coalbure.hpp
 * Project: collisions
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 24th May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to enact collision-
 * coalescence, breakup or rebound events
 * in SDM analagous to to Shima et al. 2009.
 * CoalBuRe struct satisfies PairEnactX
 * concept used in Collisions struct
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_COALBURE_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_COALBURE_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <functional>

#include "../microphysicalprocess.hpp"
#include "../superdrop.hpp"
#include "./breakup.hpp"
#include "./breakup_nfrags.hpp"
#include "./coalbure_flag.hpp"
#include "./coalescence.hpp"
#include "./collisionkinetics.hpp"
#include "./collisions.hpp"

/**
 * @brief DoCoalBuRe = DoCoalescenceBreakupRebound, i.e. enacts collision-coalescence, breakup, or
 * rebound of super-droplets.
 *
 * This class template implements the collision-coalescence, breakup, or rebound of
 * superdroplets based on specified flag values.
 *
 * @tparam NFrags Calculation for number of fragments in case of breakup.
 * @tparam Flag Flag indicating the type of action to perform: coalescence, breakup, or rebound.
 */
template <NFragments NFrags, CoalBuReFlag Flag>
struct DoCoalBuRe {
 private:
  DoCoalescence coal;   /**< Instance of DoCoalescence. */
  DoBreakup<NFrags> bu; /**< Instance of DoBreakup with specified no. of fragments calculation. */
  Flag coalbure_flag;   /**< Instance of CoalBuReFlag indicating the action to perform. */

  /**
   * @brief Calculates the value of the gamma factor in collision-coalescence.
   *
   * This function calculates the value of the gamma factor in collision-coalescence
   * based on the given probability of collision.
   *
   * _Note:_ Probability is probability of collision *NOT* collision-coalescence.
   *
   * @param xi1 Xi value of the first superdroplet.
   * @param xi2 Xi value of the second superdroplet.
   * @param prob Probability of collision.
   * @param phi Phi value.
   * @return The calculated gamma factor.
   */
  KOKKOS_FUNCTION
  uint64_t collision_gamma(const uint64_t xi1, const uint64_t xi2, const double prob,
                           const double phi) const {
    return coal.coalescence_gamma(xi1, xi2, prob, phi);
  }

  /**
   * @brief Enacts rebound, coalescence, or breakup based on the flag.
   *
   * This function enacts rebound, coalescence, or breakup based on the specified flag value:
   * If flag = 1 -> coalescence.
   * If flag = 2 -> breakup.
   * Otherwise -> rebound.
   *
   * @param gamma The gamma factor.
   * @param phi Phi value.
   * @param drop1 First superdroplet.
   * @param drop2 Second superdroplet.
   * @return True if the resulting superdroplet is null, otherwise false.
   */
  KOKKOS_FUNCTION
  bool coalesce_breakup_or_rebound(const uint64_t gamma, const double phi, Superdrop &drop1,
                                   Superdrop &drop2) const;

 public:
  /**
   * @brief Constructs a new DoCoalBuRe = DoCoalescenceBreakupRebound object.
   *
   * @param nfrags Calculation for the nmber of fragments in cases of breakup.
   * @param flag Flag indicating the action to perform: coalescence, breakup or rebound.
   */
  DoCoalBuRe(const NFrags nfrags, const Flag flag) : bu(nfrags), coalbure_flag(flag) {}

  /**
   * @brief Operator used as an adaptor such that DoCoalBuRe satisfies the PairEnactX concept
   * and so can be used as the EnactCollision function-like object in the DoCollisions struct.
   *
   * This operator calls functions to enact the collision- coalescence, breakup or rebound of
   * two super-droplets.
   *
   * @param drop1 First superdroplet.
   * @param drop2 Second superdroplet.
   * @param prob Probability of collision.
   * @param phi Phi value.
   * @return True if the resulting superdroplet is null, otherwise false.
   */
  KOKKOS_INLINE_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2, const double prob, const double phi) const;
};

/**
 * @brief Constructs a Microphysical Process for collision-coalescence, breakup, or rebound of
 * superdroplets.
 *
 * This function constructs a Microphysical Process for collision-coalescence, breakup, or rebound
 * of superdroplets with a constant timestep 'interval' and probability of collision determined by
 * 'collprob'.
 *
 * @tparam Probability Type of PairProbability.
 * @tparam NFrags Number of fragments for breakup.
 * @tparam Flag Flag indicating the action to perform.
 * @param interval Timestep interval between collision events.
 * @param int2realtime Function to convert interval to a real time [s].
 * @param collprob Probability of collisions.
 * @param nfrags Calculatino for number of fragments cases of breakup.
 * @param coalbure_flag Flag indicating the action to perform: coalescence, breakup or rebound.
 * @return A Microphysical Process enacting collision- coalescence, breakup or rebound.
 */
template <PairProbability Probability, NFragments NFrags, CoalBuReFlag Flag>
inline MicrophysicalProcess auto CoalBuRe(const unsigned int interval,
                                          const std::function<double(unsigned int)> int2realtime,
                                          const Probability collprob, const NFrags nfrags,
                                          const Flag coalbure_flag) {
  const auto DELT = double{int2realtime(interval)};

  const DoCoalBuRe<NFrags, Flag> coalbure(nfrags, coalbure_flag);
  const DoCollisions<Probability, DoCoalBuRe<NFrags, Flag>> colls(DELT, collprob, coalbure);

  return ConstTstepMicrophysics(interval, colls);
}

/**
 * @brief Operator used as an adaptor such that DoCoalBuRe satisfies the PairEnactX concept
 * and so can be used as the EnactCollision function-like object in the DoCollisions struct.
 *
 * This operator calls functions to enact the collision- coalescence, breakup or rebound of
 * two super-droplets.
 *
 * @param drop1 First superdroplet.
 * @param drop2 Second superdroplet.
 * @param prob Probability of collision.
 * @param phi Phi value.
 * @return True if the resulting superdroplet is null, otherwise false.
 */
template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool DoCoalBuRe<NFrags, Flag>::operator()(Superdrop &drop1, Superdrop &drop2,
                                                          const double prob,
                                                          const double phi) const {
  /* 1. calculate gamma factor for collision  */
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();
  const auto gamma = collision_gamma(xi1, xi2, prob, phi);

  /* 2. enact collision between pair
  of superdroplets if gamma is not zero */
  if (gamma != 0) {
    return coalesce_breakup_or_rebound(gamma, phi, drop1, drop2);
  }

  return 0;
}

/**
 * @brief Enacts rebound, coalescence, or breakup based on the flag.
 *
 * This function enacts rebound, coalescence, or breakup based on the specified flag value:
 * If flag = 1 -> coalescence.
 * If flag = 2 -> breakup.
 * Otherwise -> rebound.
 *
 * @param gamma The gamma factor.
 * @param phi Phi value.
 * @param drop1 First superdroplet.
 * @param drop2 Second superdroplet.
 * @return True if the resulting superdroplet is null, otherwise false.
 */
template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool DoCoalBuRe<NFrags, Flag>::coalesce_breakup_or_rebound(const uint64_t gamma,
                                                                           const double phi,
                                                                           Superdrop &drop1,
                                                                           Superdrop &drop2) const {
  const auto flag = coalbure_flag(phi, drop1, drop2);

  bool is_null(0);
  switch (flag) {
    case 1:  // coalescence
      is_null = coal.coalesce_superdroplet_pair(gamma, drop1, drop2);
      break;
    case 2:  // breakup
      bu.breakup_superdroplet_pair(drop1, drop2);
      break;
  }

  return is_null;
}

#endif  // LIBS_SUPERDROPS_COLLISIONS_COALBURE_HPP_
