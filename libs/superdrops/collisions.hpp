/*
 * ----- CLEO -----
 * File: collisions.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 16th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for modelling collision
 * microphysical processes in SDM
 * e.g. collision-coalescence
 */

#ifndef COLLISIONS_HPP
#define COLLISIONS_HPP

#include <concepts>

#include <Kokkos_Core.hpp>

#include "./kokkosaliases_sd.hpp"
#include "./superdrop.hpp"
#include "./state.hpp"
#include "./urbg.hpp"

template <typename P>
concept PairProbability = requires(P p,
                                   Superdrop &drop,
                                   double d)
/* Objects that are of type 'PairProbability'
take a pair of superdroplets and returns
something convertible to a double
(hopefully a probability!) */
{
  {
    p(drop, drop, d, d)
  } -> std::convertible_to<double>;
};

template <typename X>
concept PairEnactX = requires(X x,
                              Superdrop &drop,
                              double d)
/* Objects that are of type PairEnactX
takes a pair of superdrops and returns
void (it may change the properties of
the superdrops)*/
{
  {
    x(drop, drop, d, d)
  } -> std::same_as<void>;
};

template <PairProbability Probability,
          PairEnactX EnactCollision>
struct DoCollisions
/* class for method to enact collisions
between superdrops e.g. collision-coalescence */
{
private:
  const double DELT; // real time interval [s] for which probability of collision-x is calculated
  
  const Probability probability;
  /* object (has operator that) returns prob_jk, the probability
  a pair of droplets undergo some kind of collision process.
  prob_jk is analogous to prob_jk = K(drop1, drop2) delta_t/delta_vol,
  where K(drop1, drop2) := C(drop1, drop2) * |v1−v2|
  is the coalescence kernel (see Shima 2009 eqn 3). For example
  prob_jk may return the probability of collision-coalescence
  according to a particular coalescence kernel, or collision-breakup */

  const EnactCollision collision;
  /* object (has operator that) enacts a collision-X event on two
  superdroplets. For example it may enact collision-coalescence by
  of a pair of superdroplets by changing the multiplicity,
  radius and solute mass of each superdroplet in the pair
  according to Shima et al. 2009 Section 5.1.3. part (5). */

  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION
  subviewd_supers do_collisions(subviewd_supers supers,
                                const State &state,
                                URBG<DeviceType> urbg) const
  /* Superdroplet collision method adapted from collision-coalescence
  in Shima et al. 2009. This function shuffles supers to get random
  pairs of superdroplets (SDs) and then calls the collision function
  for each pair (assuming these SDs are colliding some 'VOLUME' [m^3]) */
  {
    const double VOLUME(state.get_volume() * dlc::VOL0);    // volume in which collisions occur [m^3]
    const size_t nsupers(supers.extent(0));
    const size_t nhalf(nsupers / 2); // same as floor() for positive nsupers
    const double scale_p(nsupers * (nsupers - 1.0) / (2.0 * nhalf));

    /* Randomly shuffle order of superdroplet objects
    in order to generate random pairs */
    shuffle_supers(supers, urbg);
    
    /* collide all randomly generated pairs of SDs */
    for (size_t i = 1; i < nsupers; i += 2)
    {
      collide_superdroplet_pair(supers(i - 1), supers(i),
                                urbg, scale_p, VOLUME);
    }

    // return remove_outofdomain_superdrops(span4SDsinGBx);

    return supers;
  }

  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION void
  collide_superdroplet_pair(Superdrop &dropA,
                            Superdrop &dropB,
                            URBG<DeviceType> &urbg,
                            const double scale_p,
                            const double VOLUME) const
  /* Monte Carlo Routine from Shima et al. 2009 for
  collision-coalescence generalised to any collision-X
  process for a pair of superdroplets */
  {
    /* 1. assign references to each superdrop in pair
    that will collide such that (drop1.xi) >= (drop2.xi) */
    const auto drops(assign_drops(dropA, dropB)); // {drop1, drop2}

    /* 2. calculate scaled probability of
    collision for pair of superdroplets */
    const double prob(scaled_probability(drops.first,
                                         drops.second,
                                         scale_p, VOLUME));

    /* 3. Monte Carlo Step: use random number to
    enact (or not) collision of superdroplets pair */
    const double phi(urbg.drand(0.0, 1.0)); // random number in range [0.0, 1.0]
    enact_collisionx(drops.first, drops.second, prob, phi);
  }

  KOKKOS_INLINE_FUNCTION Kokkos::pair<Superdrop &, Superdrop &>
  assign_drops(Superdrop &dropA, Superdrop &dropB)
  /* compare dropA.xi with dropB.xi and return (non-const)
  references to dropA and dropB in a pair {drop1, drop2}
  such that drop1.xi is always >= drop2.xi */
  {
    if (!(dropA.xi < dropB.xi))
    {
      return {dropA, dropB}; 
    }
    else
    {
      return {dropB, dropA};
    }
  }

  KOKKOS_INLINE_FUNCTION double
  scaled_probability(const Superdrop &drop1,
                     const Superdrop &drop2,
                     const double scale_p,
                     const double VOLUME)
  /* calculate probability of pair of superdroplets
  undergoing collision-x according to Shima et al. 2009
  ("p_alpha" in paper). Assumes drop1.xi >= drop2.xi */
  {
    const double prob_jk(probability(drop1, drop2, DELT, VOLUME));
    const double large_xi(drop1.get_xi()); // casting xi to double (!)
    
    const double prob(scale_p * large_xi * prob_jk); 

    return prob;
  } 

public:
  DoCollisions(const double DELT, Probability p, EnactCollision x)
      : DELT(DELT), probability(p), collision(x) {}

  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION
  subviewd_supers operator()(const unsigned int subt,
                            subviewd_supers supers,
                            const State &state,
                            URBG<DeviceType> urbg) const
  /* this operator is used as an "adaptor" for using
  collisions as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    do_collisions(supers, state, urbg);

    return supers;
  }
};

#endif // COLLISIONS_HPP