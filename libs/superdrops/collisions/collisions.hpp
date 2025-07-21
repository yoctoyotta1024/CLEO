/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: collisions.hpp
 * Project: collisions
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct for modelling collision microphysical processes in SDM e.g. collision-coalescence
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_COLLISIONS_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_COLLISIONS_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <concepts>
#include <random>

#include "../../cleoconstants.hpp"
#include "../kokkosaliases_sd.hpp"
#include "../sdmmonitor.hpp"
#include "../state.hpp"
#include "../superdrop.hpp"
#include "./shuffle.hpp"

namespace dlc = dimless_constants;

/**
 * @brief Concept for objects that return a probability of collision between two (real) droplets.
 *
 * Object (has operator that) returns prob_jk, the probability a pair of droplets undergo some
 * kind of collision process. Usually prob_jk = K(drop1, drop2) delta_t/delta_vol,
 * where K(drop1, drop2) := C(drop1, drop2) * |v1−v2| is the coalescence kernel
 * (see Shima 2009 eqn 3). For example a type of PairProbability may return prob_jk which is the
 * probability of collision-coalescence according to a particular coalescence kernel, or
 * collision-breakup.
 *
 * @tparam P The type representing the pair probability object.
 */
template <typename P>
concept PairProbability = requires(P p, Superdrop &drop, double d) {
  { p(drop, drop, d, d) } -> std::convertible_to<double>;
};

/**
 * @brief Concept for objects that enact a sucessful collision event between two superdroplets, e.g.
 * to model the coalscence and/or rebound and/or breakup of two superdroplets.
 *
 * Object (has operator that) enacts a collision-X event between two superdroplets. For example it
 * may enact collision-coalescence of a pair of superdroplets by changing the multiplicity,
 * radius and solute mass of each superdroplet in the pair.
 *
 * @tparam X The type representing the pair enactment object.
 */
template <typename X>
concept PairEnactX = requires(X x, Superdrop &drop, double d) {
  { x(drop, drop, d, d) } -> std::convertible_to<bool>;
};

/*
CollideSupersFunctor struct encapsulates superdroplet collisions so that parallel loop
in collide_supers function (see below) only captures necessary objects and not
other members of DoCollisions coincidentally
*/
template <PairProbability Probability, PairEnactX EnactCollision>
struct CollideSupersFunctor {
  const Probability &probability;        /**< Object for calculating collision probabilities. */
  const EnactCollision &enact_collision; /**< Enactment object for enacting collision events. */
  const GenRandomPool genpool;           /**< Kokkos thread-safe random number generator pool.*/
  const subviewd_supers supers;          /**< The randomly shuffled view of super-droplets. */
  const double scale_p;                  /**< The probability scaling factor. */
  const double DELT;   /**< time interval [s] over which probability of collision is calculated. */
  const double VOLUME; /**< The volume [m^-3]. */

  /**
   * @brief Assigns references to super-droplets in a pair based on their multiplicities.
   *
   * Compare dropA's multiplicity with dropB's, and returns (non-const) references to dropA
   * and dropB in a pair {drop1, drop2} such that drop1's multiplicity is always >= drop2's. *
   *
   * @param dropA The first super-droplet.
   * @param dropB The second super-droplet.
   * @return A pair of references to super-droplets ordered by descending xi value.
   */
  KOKKOS_INLINE_FUNCTION Kokkos::pair<Superdrop &, Superdrop &> assign_drops(
      Superdrop &dropA, Superdrop &dropB) const {
    if (!(dropA.get_xi() < dropB.get_xi())) {
      return {dropA, dropB};
    } else {
      return {dropB, dropA};
    }
  }

  /**
   * @brief Scaled probability of collision for a pair of super-droplets.
   *
   * Returns the probability of pair of super-droplets colliding according to
   * Shima et al. 2009 ("p_alpha" in paper). Function assumes drop1.xi >= drop2.xi.
   *
   * _Note:_ multiplicity, xi, of drop1 is cast to double for the calculation.
   *
   * @param drop1 The first super-droplet.
   * @param drop2 The second super-droplet.
   * @param scale_p The probability scaling factor.
   * @param VOLUME The volume [m^-3].
   * @return The scaled probability of the collision.
   */
  KOKKOS_INLINE_FUNCTION double scaled_probability(const Superdrop &drop1, const Superdrop &drop2,
                                                   const double scale_p,
                                                   const double VOLUME) const {
    const auto prob_jk = double{probability(drop1, drop2, DELT, VOLUME)};
    const auto large_xi = static_cast<double>(drop1.get_xi());  // casting to double (!)

    const auto prob = double{scale_p * large_xi * prob_jk};

    return prob;
  }

  /**
   * @brief Performs collision event for a pair of superdroplets.
   *
   * Monte Carlo Routine from Shima et al. 2009 for collision-coalescence generalised to any
   * collision-[X] process for a pair of super-droplets.
   *
   * @param dropA The first superdroplet.
   * @param dropB The second superdroplet.
   * @param scale_p The probability scaling factor.
   * @param VOLUME The volume [m^-3].
   * @return True if the collision event results in null superdrops with xi=0), otherwise false.
   */
  KOKKOS_INLINE_FUNCTION void collide_superdroplet_pair(Superdrop &dropA, Superdrop &dropB,
                                                        const double scale_p,
                                                        const double VOLUME) const {
    /* 1. assign references to each superdrop in pair that will collide
    such that (drop1.xi) >= (drop2.xi) */
    const auto drops = assign_drops(dropA, dropB);  // {drop1, drop2}

    /* 2. calculate scaled probability of collision for pair of superdroplets */
    const auto prob = scaled_probability(drops.first, drops.second, scale_p, VOLUME);

    /* 3. Monte Carlo Step: use random number to enact (or not) collision of superdroplets pair */
    URBG<ExecSpace> urbg{genpool.get_state()};  // thread safe random number generator
    const auto phi = urbg.drand(0.0, 1.0);      // random number in range [0.0, 1.0]
    genpool.free_state(urbg.gen);

    enact_collision(drops.first, drops.second, prob, phi);
  }

  /*
   * operator for functor with parallel (TeamThreadRangePolicy) loop over superdroplet pairs
   * in supers view in order to call collide_superdroplet_pair
   */
  KOKKOS_INLINE_FUNCTION void operator()(const size_t jj) const {
    const auto kk = size_t{jj * 2};
    collide_superdroplet_pair(supers(kk), supers(kk + 1), scale_p, VOLUME);
  }
};

/**
 * @struct DoCollisions
 * @brief Implements microphysical processes for collisions between superdroplets.
 * @tparam Probability The type representing the pair probability object.
 * @tparam EnactCollision The type representing the pair enactment object.
 */
template <PairProbability Probability, PairEnactX EnactCollision>
struct DoCollisions {
 private:
  double DELT; /**< time interval [s] over which probability of collision is calculated. */
  Probability probability; /**< Probability object for calculating collision probabilities. */
  EnactCollision enact_collision; /**< Enactment object for enacting collision events. */
  GenRandomPool genpool;          /**< Kokkos thread-safe random number generator pool.*/

  /**
   * @brief Performs collisions between super-droplets in supers view.
   *
   * Enacts collisions for pairs of super-droplets in supers view adapted from collision-coalescence
   * of Shima et al. 2009 to generalise to allow for other types of collision-[X] events.
   *
   * Function uses Kokkos nested parallelism for paralelism over supers inside parallelised loop
   * for member 'teamMember'.
   * In serial Kokkos::parallel_for([...]) is equivalent to loop:
   * for (size_t jj(0); jj < npairs; ++jj) {[...]}.
   *
   * _NOTE:_ function assumes supers is already randomly shuffled and these superdrops are
   * colliding some 'VOLUME' [m^3]).
   *
   * @param team_member The Kokkos team member.
   * @param supers The randomly shuffled view of super-droplets.
   * @param volume The volume in which to calculate the probability of collisions.
   * @return The number of null (xi=0) superdrops.
   */
  KOKKOS_INLINE_FUNCTION void collide_supers(const TeamMember &team_member, subviewd_supers supers,
                                             const double volume) const {
    const auto nsupers = static_cast<size_t>(supers.extent(0));
    const auto npairs = size_t{nsupers / 2};  // no. pairs of superdrops (=floor() for nsupers > 0)
    const auto scale_p = double{nsupers * (nsupers - 1.0) / (2.0 * npairs)};
    const auto VOLUME = double{volume * dlc::VOL0};  // volume in which collisions occur [m^3]

    const auto functor =
        CollideSupersFunctor{probability, enact_collision, genpool, supers, scale_p, DELT, VOLUME};
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, npairs), functor);
    team_member.team_barrier();  // synchronise threads
  }

  /**
   * @brief Executes collision events for pairs of superdroplets.
   *
   * Superdroplet collision algorithm adapted from collision-coalescence in Shima et al. 2009.
   * This function shuffles supers to get random pairs of superdroplets (SDs) and then calls the
   * collision function for each pair assuming these superdrops are colliding some 'VOLUME' [m^3].
   * Function is designed to be called inside a parallelised loop for member 'teamMember'.
   *
   * @param team_member The Kokkos team member.
   * @param supers The view of super-droplets.
   * @param volume The volume in which to calculate the probability of collisions.
   * @return The updated superdroplets.
   */
  KOKKOS_INLINE_FUNCTION void do_collisions(const TeamMember &team_member, subviewd_supers supers,
                                            const double volume) const {
    /* Randomly shuffle order of superdroplet objects
    in supers in order to generate random pairs */
    supers = shuffle_supers(team_member, supers, genpool);

    /* collide all randomly generated pairs of SDs */
    collide_supers(team_member, supers, volume);
  }

 public:
  /**
   * @brief Constructs a DoCollisions object.
   *
   * _Note:_ If DoCollisions used at the MicrophysicsFunction type for a ConstTstepMicrophysics
   * instance, the interval between calls of DoCollisions operator() in model timesteps must be
   * concordant with DELT [s].
   *
   * @param DELT Time interval [s] over which probability of collision is calculated.
   * @param p The probability object for calculating the probability of a collision.
   * @param x The enactment object for enacting collision events.
   */
  DoCollisions(const double DELT, Probability p, EnactCollision x)
      : DELT(DELT), probability(p), enact_collision(x), genpool(std::random_device {}()) {}

  /**
   * @brief Operator used as an "adaptor" for using collisions as the MicrophysicsFunction type for
   * a ConstTstepMicrophysics instance (*hint* which itself satsifies the MicrophysicalProcess
   * concept).
   *
   * i.e. Operator allows DoCollisions to be used as the function in a microphysical process with a
   * constant timestep between events. _Note:_ If object used in this way, the interval between
   * calls to this function (i.e. between collision events) in model timesteps should be concordant
   * witih DELT of the instance.
   *
   * @param team_member The Kokkos team member.
   * @param subt The sub-time step.
   * @param supers The superdroplets.
   * @param state The state.
   * @param mo Monitor of SDM processes.
   * @return The updated superdroplets.
   */
  KOKKOS_INLINE_FUNCTION void operator()(const TeamMember &team_member, const unsigned int subt,
                                         subviewd_supers supers, const State &state,
                                         const SDMMonitor auto mo) const {
    do_collisions(team_member, supers, state.get_volume());
  }
};

#endif  // LIBS_SUPERDROPS_COLLISIONS_COLLISIONS_HPP_
