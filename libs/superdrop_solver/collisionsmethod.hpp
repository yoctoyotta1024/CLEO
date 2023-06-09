// Author: Clara Bayley
// File: collisionsmethod.hpp
/* Header file for class that controls
collision events in superdroplet model */

#ifndef COLLISIONSMETHOD_HPP
#define COLLISIONSMETHOD_HPP

#include <concepts>
#include <random>
#include <algorithm>
#include <string>
#include <iostream>
#include <stdexcept>
#include <span>
#include <cmath>
#include <functional>
#include <concepts>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./thermostate.hpp"
#include "./coalescencekernel.hpp"
#include "./sdmprocess.hpp"
#include "./randomgen.hpp"

namespace dlc = dimless_constants;

template <typename P>
concept PairProbability = requires(P p, const Superdrop &d1,
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

template <PairProbability PairCoalescenceProbability>
class CollisionsMethod
/* class for method to enact collisions between
superdrops during collision events in SDM */
{
private:
  const double DELT; // time interval [s] for which probability of coalescence is calculated

  const PairCoalescenceProbability pair_coalesce_probability;
  /* object (has operator () that) returns probability a pair of
  droplets coalesces according to a particular coalescence kernel.
  Equation is: prob_jk = K(drop1, drop2) delta_t/delta_vol
  where K(drop1, drop2) := C(drop1, drop2) * |v1−v2|,
  is coalescence kernel (see Shima 2009 eqn 3) */

  template <class DeviceType>
  void collide_superdroplets(std::span<SuperdropWithGbxindex> span4SDsinGBx,
                             URBG<DeviceType> &urbg, const double VOLUME) const
  /* Superdroplet collision-coalescence method according to Shima et al. 2009.
  For some 'VOLUME' [m^3] in which the collisions occur, this function
  determines whether of not coalescence occurs from Monte-Carlo collisions of
  random pairs of SDs. If coalescene occurs between two superdrops, it then
  also changes the multiplicity, radius and solute mass of the superdroplets
  that coalesce. In this implementation, superdroplets are stored contiguously
  in memory in a vector composed of SuperdropWithGbxindex instances. The span
  points to the section of that vector containing the superdroplets involed
  in the collision event. */
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
      collide_superdroplet_pair(urbg, span4SDsinGBx[i-1].superdrop,
                                span4SDsinGBx[i].superdrop, scale_p,
                                VOLUME);
    }
 
  }

  template <class DeviceType>
  void collide_superdroplet_pair(URBG<DeviceType> &urbg, Superdrop &dropA,
                                 Superdrop &dropB, const double scale_p,
                                 const double VOLUME) const
  /* Monte Carlo Routine according to Shima et al. 2009
  for collision-coalescence of a pair of superdroplets.
  (Requires CollisionMethod instance has an object
  satisfying the concept of a PairCoalescenceProbability) */
  {
    /* 1. assign references to each superdrop in pair
    that will collide such that (drop1.eps) >= (drop2.eps) */
    Superdrop &drop1 = assign_superdroplet(dropA, dropB, "drop1");
    Superdrop &drop2 = assign_superdroplet(dropA, dropB, "drop2");

    /* make copies of eps1 and eps2 for ease of use */
    const unsigned long long eps1 = drop1.eps;
    const unsigned long long eps2 = drop2.eps;

    /* 2. determine scaled probability of pair coalescence
    according to Shima et al. 2009 ("p_alpha" in paper) */
    const double prob_jk = pair_coalesce_probability(drop1, drop2, DELT, VOLUME);
    const double prob = scale_p * std::max(eps1, eps2) * prob_jk;

    /* 3. Monte Carlo step: randomly determine coalescence gamma factor */
    const unsigned long long gamma = monte_carlo_gamma(urbg, prob, eps1, eps2);

    /* 4. coalesce particles if gamma != 0 */
    if (gamma != 0)
    {
      coalesce_superdroplet_pair(drop1, drop2, gamma);
    }
  }

  Superdrop &assign_superdroplet(Superdrop &dropA, Superdrop &dropB,
                                 const std::string whichdrop) const
  /* compare dropA.eps with dropB.eps and return
  either drop1 or drop2 such that drop1.eps is always >
  drop2.eps */
  {
    if (dropA.eps > dropB.eps)
    {
      if (whichdrop == "drop1")
      {
        return dropA;
      }
      else
      {
        return dropB;
      }
    }

    else
    {
      if (whichdrop == "drop1")
      {
        return dropB;
      }
      else
      {
        return dropA;
      }
    }
  }

  template <class DeviceType>
  unsigned long long monte_carlo_gamma(URBG<DeviceType> &urbg, const double prob,
                                       const unsigned long long eps1,
                                       const unsigned long long eps2) const
  /* calculates value of gamma factor in
  Monte Carlo collision-coalescence process
  according to Shima et al. 2009 */
  {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    const double phi = dis(urbg); // random number phi in range [0,1]

    unsigned long long gamma = 0;
    if (phi < (prob - floor(prob)))
    {
      gamma = floor(prob) + 1;
    }
    else if (phi >= (prob - floor(prob)))
    {
      gamma = floor(prob);
    }

    const unsigned long long maxgamma(eps1 / eps2); // same as floor() for positive ints

    return std::min(gamma, maxgamma);
  }

  void coalesce_superdroplet_pair(Superdrop &drop1, Superdrop &drop2,
                                  const unsigned long long gamma) const
  /* coalesce pair of superdroplets by changing multiplicity,
  radius and solute mass of each superdroplet in pair
  according to Shima et al. 2009 Section 5.1.3. part (5) */
  {
    if (drop1.eps - gamma * drop2.eps == 0)
    {
      twin_superdroplet_coalescence(drop1, drop2, gamma);
    }

    else if (drop1.eps - gamma * drop2.eps > 0)
    {
      different_superdroplet_coalescence(drop1, drop2, gamma);
    }

    else
    {
      std::string errormsg = "something undefined occured during colllision-coalescence" +
                             std::to_string(drop1.eps) + " < " +
                             std::to_string(gamma * (drop2.eps)) + " ?";
      throw std::invalid_argument(errormsg);
    }
  }

  void twin_superdroplet_coalescence(Superdrop &drop1,
                                     Superdrop &drop2,
                                     const unsigned long long gamma) const
  /* if eps1 = gamma*eps2 coalescence makes twin SDs
  with same eps, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */
  { 
    const unsigned long long new_eps = drop2.eps / 2;
    const double new_m_sol = drop2.m_sol + gamma * drop1.m_sol;
    const double new_rcubed = pow(drop2.radius, 3.0) + gamma * pow(drop1.radius, 3.0);
    const double new_r = pow(new_rcubed, (1.0 / 3.0));

    drop1.eps = new_eps;
    drop2.eps = drop2.eps - new_eps;

    drop1.radius = new_r;
    drop2.radius = new_r;

    drop1.m_sol = new_m_sol;
    drop2.m_sol = new_m_sol;
  }

  void different_superdroplet_coalescence(Superdrop &drop1,
                                          Superdrop &drop2,
                                          const unsigned long long gamma) const
  /* if eps1 > gamma*eps2 coalescence grows drop2 radius and mass
  via decreasing multiplicity of drop1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */
  {
    drop1.eps = drop1.eps - gamma * drop2.eps;

    const double new_rcubed = pow(drop2.radius, 3.0) + gamma * pow(drop1.radius, 3.0);
    drop2.radius = pow(new_rcubed, (1.0 / 3.0));
    drop2.m_sol = drop2.m_sol + gamma * drop1.m_sol;
  }

public:
  CollisionsMethod(const double DELT, PairCoalescenceProbability p)
      : DELT(DELT), pair_coalesce_probability(p) {}

  template <class DeviceType>
  inline void operator()(const int currenttimestep,
                  std::span<SuperdropWithGbxindex> span4SDsinGBx,
                  ThermoState &state,
                  URBG<DeviceType> &urbg) const
  /* this operator is used as an "adaptor" for using a run_step
  function in order to call collide_superdroplets. (*hint* run_step
  is usually found within a type that satisfies the SdmProcess concept) */
  {
    const double VOLUME = state.get_volume() * pow(dlc::COORD0, 3.0); // volume in which collisions occur [m^3]
    collide_superdroplets(span4SDsinGBx, urbg, VOLUME);
  }
};

template <PairProbability PairCoalescenceProbability>
SdmProcess auto CollisionsProcess(const int interval,
                                  const std::function<double(int)> int2time,
                                  const PairCoalescenceProbability p)
{
  const double realtstep = int2time(interval);
  return ConstTstepProcess{interval, CollisionsMethod(realtstep, p)};
}

#endif // COLLISIONSMETHOD_HPP