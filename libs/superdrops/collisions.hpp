/*
 * ----- CLEO -----
 * File: collisions.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
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
                              double p)
/* Objects that are of type PairEnactX
takes a pair of superdrops and returns
void (it may change the properties of
the superdrops)*/
{
  {
    x(drop, drop, p, p)
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

  const EnactCollision enact_collision;
  /* object (has operator that) enacts a collision-X event on two
  superdroplets. For example it may enact collision-coalescence by
  of a pair of superdroplets by changing the multiplicity,
  radius and solute mass of each superdroplet in the pair
  according to Shima et al. 2009 Section 5.1.3. part (5). */

  KOKKOS_INLINE_FUNCTION
  subviewd_supers do_collisions(const unsigned int subt,
                                const subviewd_supers supers) const
  {
    // TODO
    // const double VOLUME = state.get_volume() * pow(dlc::COORD0, 3.0); // volume in which collisions occur [m^3]
    // collide_superdroplets(span4SDsinGBx, urbg, VOLUME);

    // return remove_outofdomain_superdrops(span4SDsinGBx);

    return supers;
  }

public:
  DoCollisions(const double DELT, Probability p, EnactCollision x)
      : DELT(DELT), probability(p), collision(x) {}

  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION
  subviewd_supers operator()(const unsigned int subt,
                            subviewd_supers supers,
                            State &state,
                            URBG<DeviceType> &urbg) const
  /* this operator is used as an "adaptor" for using
  collisions as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    do_collisions(subt, supers);

    return supers;
  }
};

#endif // COLLISIONS_HPP