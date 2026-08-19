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
#include <cstdint>
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
   * @brief Calculates the value of the gamma factor in Monte Carlo collision.
   *
   * This function calculates the value of the gamma factor for collisions
   * based on the given probability of collisions. The calculation is as described for
   * collision-coalescence in Shima et al. 2009 but applied to collisions (which may result in
   * coalescence, rebound or breakup) not just collision-coalescence.
   *
   * _Note:_ Probability is probability of collision *NOT* collision-coalescence.
   *
   * @param xi1 The multiplicity of the first super-droplet.
   * @param xi2 The multiplicity of the second super-droplet.
   * @param prob The probability of collision.
   * @param phi Random number in the range [0.0, 1.0].
   * @return The calculated value of the collision gamma factor.
   */
  KOKKOS_FUNCTION
  uint64_t collision_gamma(const uint64_t xi1, const uint64_t xi2, const double prob,
                           const double phi) const {
    uint64_t gamma = Kokkos::floor(prob);  // if phi >= (prob - floor(prob))
    if (phi < (prob - gamma)) {
      ++gamma;
    }

    const auto maxgamma = xi1 / xi2;  // same as floor() for positive ints

    return Kokkos::fmin(gamma, maxgamma);
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
   * @return Resulting total number of null (xi=0) superdroplets from collision (i.e. 0, 1 or 2)
   */
  KOKKOS_FUNCTION
  size_t coalesce_breakup_or_rebound(const uint64_t gamma, const double phi, Superdrop& drop1,
                                     Superdrop& drop2) const;

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
   * *NOTE:* phi_out is used here for (second) Monte Carlo step, to determine if outcome of a
   * collision is coalescene, rebound or breakup. (re-scaling could instead be implemented to use
   * phi_coll for both collision and outcome of collision but this is not done
   * here).
   *
   * @param drop1 First superdroplet.
   * @param drop2 Second superdroplet.
   * @param prob Probability of collision.
   * @param phi_coll Random number in the range [0.0, 1.0] for collision.
   * @param phi_out Random number in the range [0.0, 1.0] for outcome of collision as breakup,
   * rebound or coalescence.
   * @return Resulting total number of null (xi=0) superdroplets from collision (i.e. 0, 1 or 2)
   */
  KOKKOS_INLINE_FUNCTION
  size_t operator()(Superdrop& drop1, Superdrop& drop2, const double prob, const double phi_coll,
                    const double phi_out) const;
};

/**
 * @brief Constructs a Microphysical Process for collision-coalescence, breakup, or rebound of
 * superdroplets.
 *
 * This function constructs a Microphysical Process for collision-coalescence, breakup, or rebound
 * of superdroplets with a constant timestep 'interval' and probability of collision determined by
 * 'collprob' and a random seed for the random number generator.
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
  const MicrophysicsFunc auto colls =
      DoCollisions<Probability, DoCoalBuRe<NFrags, Flag>>(DELT, collprob, coalbure);

  return ConstTstepMicrophysics(interval, colls);
}

/**
 * same as CoalBuRe above but with fixed seed
 */
template <PairProbability Probability, NFragments NFrags, CoalBuReFlag Flag>
inline MicrophysicalProcess auto CoalBuRe(const unsigned int interval,
                                          const std::function<double(unsigned int)> int2realtime,
                                          const Probability collprob, const NFrags nfrags,
                                          const Flag coalbure_flag, const uint64_t seed) {
  const auto DELT = double{int2realtime(interval)};

  const DoCoalBuRe<NFrags, Flag> coalbure(nfrags, coalbure_flag);
  const MicrophysicsFunc auto colls =
      DoCollisions<Probability, DoCoalBuRe<NFrags, Flag>>(DELT, collprob, coalbure, seed);

  return ConstTstepMicrophysics(interval, colls);
}

/**
 * @brief Operator used as an adaptor such that DoCoalBuRe satisfies the PairEnactX concept
 * and so can be used as the EnactCollision function-like object in the DoCollisions struct.
 *
 * This operator calls functions to enact the collision- coalescence, breakup or rebound of
 * two super-droplets.
 *
 * *NOTE:* phi_out is used here for (second) Monte Carlo step, to determine if outcome of a
 * collision is coalescene, rebound or breakup. (re-scaling could instead be implemented to use
 * phi_coll for both collision and outcome of collision but this is not done
 * here).
 *
 * @param drop1 First superdroplet.
 * @param drop2 Second superdroplet.
 * @param prob Probability of collision.
 * @param phi_coll Random number in the range [0.0, 1.0] for collision.
 * @param phi_out Random number in the range [0.0, 1.0] for outcome of collision as breakup,
 * rebound or coalescence.
 * @return Resulting total number of null (xi=0) superdroplets from collision (i.e. 0, 1 or 2)
 */
template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION size_t DoCoalBuRe<NFrags, Flag>::operator()(Superdrop& drop1, Superdrop& drop2,
                                                            const double prob,
                                                            const double phi_coll,
                                                            const double phi_out) const {
  /* 1. calculate gamma factor for collision  */
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();
  const auto gamma = collision_gamma(xi1, xi2, prob, phi_coll);

  /* 2. enact collision between pair
  of superdroplets if gamma is not zero */
  if (gamma != 0) {
    return coalesce_breakup_or_rebound(gamma, phi_out, drop1, drop2);
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
 * @return Resulting total number of null (xi=0) superdroplets from collision (i.e. 0, 1 or 2)
 */
template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION size_t DoCoalBuRe<NFrags, Flag>::coalesce_breakup_or_rebound(
    const uint64_t gamma, const double phi, Superdrop& drop1, Superdrop& drop2) const {
  const auto flag = coalbure_flag(phi, drop1, drop2);

  size_t null_sds = 0;
  switch (flag) {
    case 1:  // coalescence
      null_sds = coal.coalesce_superdroplet_pair(gamma, drop1, drop2);
      break;
    case 2:  // breakup
      null_sds = bu.breakup_superdroplet_pair(drop1, drop2);
      break;
  }

  return null_sds;
}

#endif  // LIBS_SUPERDROPS_COLLISIONS_COALBURE_HPP_
