/*
 * ----- CLEO -----
 * File: collisions.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 21st December 2023
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
#include <Kokkos_Random.hpp>

#include "../cleoconstants.hpp"
#include "./kokkosaliases_sd.hpp"
#include "./nullsuperdrops.hpp"
#include "./superdrop.hpp"
#include "./state.hpp"
#include "./urbg.hpp"

namespace dlc = dimless_constants;

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
  } -> std::convertible_to<bool>;
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

  const EnactCollision enact_collision;
  /* object (has operator that) enacts a collision-X event on two
  superdroplets. For example it may enact collision-coalescence by
  of a pair of superdroplets by changing the multiplicity,
  radius and solute mass of each superdroplet in the pair
  according to Shima et al. 2009 Section 5.1.3. part (5). */

  KOKKOS_INLINE_FUNCTION double
  scaled_probability(const Superdrop &drop1,
                     const Superdrop &drop2,
                     const double scale_p,
                     const double VOLUME) const
  /* calculate probability of pair of superdroplets
  undergoing collision-x according to Shima et al. 2009
  ("p_alpha" in paper). Assumes drop1.xi >= drop2.xi */
  {
    const auto prob_jk = double{probability(drop1, drop2, DELT, VOLUME)};
    const auto large_xi = (double)drop1.get_xi(); // casting xi to double (!)

    const auto prob = double{scale_p * large_xi * prob_jk};

    return prob;
  }

  KOKKOS_INLINE_FUNCTION Kokkos::pair<Superdrop &, Superdrop &>
  assign_drops(Superdrop &dropA, Superdrop &dropB) const
  /* compare dropA.xi with dropB.xi and return (non-const)
  references to dropA and dropB in a pair {drop1, drop2}
  such that drop1.xi is always >= drop2.xi */
  {
    if (!(dropA.get_xi() < dropB.get_xi()))
    {
      return {dropA, dropB};
    }
    else
    {
      return {dropB, dropA};
    }
  }

  KOKKOS_INLINE_FUNCTION bool
  collide_superdroplet_pair(Superdrop &dropA,
                            Superdrop &dropB,
                            GenRandomPool genpool,
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
    const auto prob = scaled_probability(drops.first,
                                         drops.second,
                                         scale_p, VOLUME);

    /* 3. Monte Carlo Step: use random number to
    enact (or not) collision of superdroplets pair */
    URBG<ExecSpace> urbg{genpool.get_state()}; // thread safe random number generator
    const auto phi = urbg.drand(0.0, 1.0); // random number in range [0.0, 1.0]
    genpool.free_state(urbg.gen);

    const auto isnull = enact_collision(drops.first,
                                        drops.second,
                                        prob, phi);

    return isnull;
  }

  KOKKOS_INLINE_FUNCTION size_t
  collide_supers(const TeamMember &team_member,
                 subviewd_supers supers,
                 const double volume, 
                 GenRandomPool genpool) const
  /* Enacts collisions for pairs of superdroplets in supers
  like for collision-coalescence in Shima et al. 2009.
  Assumes supers is already randomly shuffled and these
  superdrops are colliding some 'VOLUME' [m^3]). Function
  uses Kokkos nested parallelism for paralelism over supers
  inside parallelised loop for member 'teamMember'. In serial
  Kokkos::parallel_reduce([...]) is equivalent to summing nnull
  over for loop: for (size_t jj(0); jj < npairs; ++jj) {[...]} 
  when in serial */
  {
    const auto nsupers = (size_t)supers.extent(0);
    const auto npairs = size_t{nsupers / 2}; // no. pairs of superdroplets (same as floor() for positive nsupers)
    const auto scale_p = double{nsupers * (nsupers - 1.0) / (2.0 * npairs)};
    const auto VOLUME = double{volume * dlc::VOL0}; // volume in which collisions occur [m^3]

    auto totnnull = size_t{0}; // number of null superdrops
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team_member, npairs),
        [&, this](const size_t jj, size_t &nnull)
        {
          const auto kk = size_t{jj * 2};
          const auto isnull = collide_superdroplet_pair(supers(kk),
                                                        supers(kk + 1),
                                                        genpool,
                                                        scale_p,
                                                        VOLUME);
          nnull += (size_t)isnull;
        },
        totnnull);

    return totnnull;
  }

  KOKKOS_INLINE_FUNCTION subviewd_supers
  do_collisions(const TeamMember &team_member,
                subviewd_supers supers,
                const double volume,
                GenRandomPool genpool) const
  /* Superdroplet collision method adapted from collision-coalescence
  in Shima et al. 2009. This function shuffles supers to get
  random pairs of superdroplets (SDs) and then calls the
  collision function for each pair (assuming these superdrops
  are colliding some 'VOLUME' [m^3]). Function is called inside
  a parallelised loop for member 'teamMember'. */
  {
    /* Randomly shuffle order of superdroplet objects
    in supers in order to generate random pairs */
    supers = one_shuffle_supers(team_member, supers, genpool);

    /* collide all randomly generated pairs of SDs */
    size_t nnull(collide_supers(team_member, supers, volume, genpool)); // number of null superdrops
  
    return is_null_supers(supers, nnull);
  }

public:
  DoCollisions(const double DELT, Probability p, EnactCollision x)
      : DELT(DELT), probability(p), enact_collision(x) {}

  KOKKOS_INLINE_FUNCTION subviewd_supers
      operator()(const TeamMember &team_member,
                 const unsigned int subt,
                 subviewd_supers supers,
                 const State &state,
                 GenRandomPool genpool) const
  /* this operator is used as an "adaptor" for using
  collisions as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    return do_collisions(team_member, supers,
                         state.get_volume(), genpool);
  }
};

#endif // COLLISIONS_HPP