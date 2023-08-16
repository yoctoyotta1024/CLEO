// Author: Clara Bayley
// File: collisionx.hpp
/* Header file for class that controls
collision-[X] events in superdroplet
model, e.g. collision-coalescences or
collision-breakups */

#ifndef COLLISIONX_HPP
#define COLLISIONX_HPP

#include <algorithm>
#include <cmath>
#include <concepts>
#include <iostream>
#include <functional>
#include <random>
#include <span>
#include <string>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./thermostate.hpp"
#include "./sdmprocess.hpp"
#include "./randomgen.hpp"

namespace dlc = dimless_constants;

template <typename P>
concept SDPairProbability = requires(P p,
                                     const Superdrop &d1,
                                     const Superdrop &d2,
                                     const double t,
                                     const double v)
/* Objects that are of type 'PairProbability'
take a pair of superdroplets and returns
something convertible to a double
(hopefully a probability!) */
{
  {
    p(d1, d2, t, v)
  } -> std::convertible_to<double>;
};

template <typename X>
concept SDinGBxPairEnactX = requires(X x,
                                SuperdropWithGbxindex &SDinGBx1,
                                SuperdropWithGbxindex &SDinGBx2,
                                const double p,
                                const double i)
/* Objects that are of type SDinGBxPairEnactX 
takes a pair of superdrops and returns
void (it may change the properties of
the superdrops)*/
{
  {
    x(SDinGBx1, SDinGBx2, p, i)
    x(SDinGBx1, SDinGBx2, p, i)
  } -> std::same_as<void>;
};

template <SDPairProbability CollisionXProbability,
          SDinGBxPairEnactX CollisionXEnactment>
          SDinGBxPairEnactX CollisionXEnactment>
class CollisionX
/* class for method to enact collisions between
superdrops during collision events in SDM */
{
private:
  const double DELT; // time interval [s] for which probability of collision-x is calculated

  const CollisionXProbability collisionx_probability;
  /* object (has operator that) returns prob_jk, the probability
  a pair of droplets undergo some kind of collision process.
  prob_jk is analogous to prob_jk = K(drop1, drop2) delta_t/delta_vol,
  where K(drop1, drop2) := C(drop1, drop2) * |v1−v2|
  is the coalescence kernel (see Shima 2009 eqn 3). For example
  prob_jk may return the probability of collision-coalescence
  according to a particular coalescence kernel, or collision-breakup */

  const CollisionXEnactment enact_collisionx;
  /* object (has operator that) enacts a collision-X event on two
  superdroplets. For example it may enact collision-coalescence by
  of a pair of superdroplets by changing the multiplicity,
  radius and solute mass of each superdroplet in the pair
  according to Shima et al. 2009 Section 5.1.3. part (5). */

  template <class DeviceType>
  void collide_superdroplets(std::span<SuperdropWithGbxindex> span4SDsinGBx,
                             URBG<DeviceType> &urbg, const double VOLUME) const
  /* Superdroplet collision method adapted from collision-coalescence in
  Shima et al. 2009. This function determines the random pairs of
  superdroplets (SDs) from the span4SDsinGBx and calls the collision function
  for each pair (assuming these SDs are colliding some 'VOLUME' [m^3]) */
  {
    const size_t nsupers(span4SDsinGBx.size());
    const size_t nhalf(nsupers / 2); // same as floor() positive nsupers
    const double scale_p(nsupers * (nsupers - 1.0) / (2.0 * nhalf));

    /* Randomly shuffle order of superdroplet objects
    in order to generate random pairs */
    std::shuffle(span4SDsinGBx.begin(), span4SDsinGBx.end(), urbg);

    /* collide all randomly generated pairs of SDs */
    for (size_t i = 1; i < nsupers; i += 2)
    {
      collide_superdroplet_pair(urbg, span4SDsinGBx[i - 1],
                                span4SDsinGBx[i], scale_p,
                                VOLUME);
    }
  }

  template <class DeviceType>
  void collide_superdroplet_pair(URBG<DeviceType> &urbg,
                                 SuperdropWithGbxindex &SDinGBxA,
                                 SuperdropWithGbxindex &SDinGBxB,
                                 const double scale_p,
  void collide_superdroplet_pair(URBG<DeviceType> &urbg,
                                 SuperdropWithGbxindex &SDinGBxA,
                                 SuperdropWithGbxindex &SDinGBxB,
                                 const double scale_p,
                                 const double VOLUME) const
  /* Monte Carlo Routine from Shima et al. 2009 for
  collision-coalescence generalised to any collision-X process
  for a pair of superdroplets given the collisionx_probability
  and collisionx_superdroplet_pair */
  {
    /* 1. assign references to each superdrop in pair
    that will collide such that (sd1.eps) >= (sd2.eps) */
    auto [SDinGBx1, SDinGBx2] = assign_superdroplet(SDinGBxA, SDinGBxB);

    /* 2. calculate scaled probability of pair collision-x
    according to Shima et al. 2009 ("p_alpha" in paper) */
    const double prob_jk = collisionx_probability(SDinGBx1.superdrop,
                                                  SDinGBx2.superdrop,
                                                  DELT, VOLUME);
    const double prob = scale_p * SDinGBx1.superdrop.eps * prob_jk; // note: eps1 >= eps2 here

    /* 3. Monte Carlo Step: use random number to enact (or not)
    collision-x on pair of superdroplets */
    std::uniform_real_distribution<> dis(0.0, 1.0);
    const double phi = dis(urbg); // phi is random number in range [0,1]
    enact_collisionx(SDinGBx1, SDinGBx2, prob, phi);
    enact_collisionx(SDinGBx1, SDinGBx2, prob, phi);
  }

  std::pair<SuperdropWithGbxindex &, SuperdropWithGbxindex &>
  assign_superdroplet(const SuperdropWithGbxindex &SDinGBxA,
                      const SuperdropWithGbxindex &SDinGBxB) const
  /* compare superdropA.eps with superdropB.eps and return (non-const)
  references to SDinGBxA and SDinGBxB in a pair {SDinGBx1, SDinGBx2}
  such that superdrop1.eps is always > superdrop2.eps */
  {
    auto compare = [](const SuperdropWithGbxindex &SDinGBxA,
                      const SuperdropWithGbxindex &SDinGBxB)
    auto compare = [](const SuperdropWithGbxindex &SDinGBxA,
                      const SuperdropWithGbxindex &SDinGBxB)
    {
      return SDinGBxA.superdrop.eps < SDinGBxB.superdrop.eps; // returns true if epsA < epsB
      return SDinGBxA.superdrop.eps < SDinGBxB.superdrop.eps; // returns true if epsA < epsB
    };

    auto [SDinGBx2, SDinGBx1] = std::minmax(SDinGBxA, SDinGBxB, compare); // drop2.eps =< drop1.eps
    auto [SDinGBx2, SDinGBx1] = std::minmax(SDinGBxA, SDinGBxB, compare); // drop2.eps =< drop1.eps

    return {const_cast<SuperdropWithGbxindex &>(SDinGBx1),
            const_cast<SuperdropWithGbxindex &>(SDinGBx2)};
    return {const_cast<SuperdropWithGbxindex &>(SDinGBx1),
            const_cast<SuperdropWithGbxindex &>(SDinGBx2)};
  }

public:
  CollisionX(const double DELT,
             CollisionXProbability p,
             CollisionXEnactment x)
      : DELT(DELT),
        collisionx_probability(p),
        enact_collisionx(x) {}

  template <class DeviceType>
  auto operator()(const int currenttimestep, const unsigned int gbxindex,
             std::span<SuperdropWithGbxindex> span4SDsinGBx,
             ThermoState &state,
             URBG<DeviceType> &urbg) const
  /* this operator is used as an "adaptor" for using a run_step
  function in order to call collide_superdroplets. (*hint* run_step
  is usually found within a type that satisfies the SdmProcess concept) */
  {
    const double VOLUME = state.get_volume() * pow(dlc::COORD0, 3.0); // volume in which collisions occur [m^3]
    collide_superdroplets(span4SDsinGBx, urbg, VOLUME);

    return remove_outofdomain_superdrops(span4SDsinGBx);
  }
};

#endif // COLLISIONX_HPP