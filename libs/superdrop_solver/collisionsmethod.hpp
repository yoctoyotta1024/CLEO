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

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./thermostate.hpp"
#include "./coalescencekernel.hpp"

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

  PairCoalescenceProbability pair_coalesce_probability;
  /* object (has operator () that) returns probability a pair of
  droplets coalesces according to a particular coalescence kernel.
  Equation is: prob_jk = K(drop1, drop2) delta_t/delta_vol
  where K(drop1, drop2) := C(drop1, drop2) * |v1−v2|,
  is coalescence kernel (see Shima 2009 eqn 3) */

  void collide_superdroplets(std::span<SuperdropWithGridbox> span4SDsinGBx,
                             std::mt19937 &gen,
                             const double VOLUME) const
  /* Superdroplet collision-coalescence method according to Shima et al. 2009.
  For some 'VOLUME' [m^3] in which the collisions occur, this function
  determines whether of not coalescence occurs from Monte-Carlo collisions of
  random pairs of SDs. If coalescene occurs between two superdrops, it then
  also changes the multiplicity, radius and solute mass of the superdroplets
  that coalesce. In this implementation, superdroplets are stored contiguously
  in memory in a vector composed of SuperdropWithGridbox instances. The span
  points to the section of that vector containing the superdroplets involed
  in the collision event. */
  {
    const int nsupers = span4SDsinGBx.size();
    const int nhalf = floor(nsupers / 2.0);
    const int scale_p = nsupers * (nsupers - 1.0) / (2.0 * nhalf);

    /* Randomly shuffle order of superdroplet objects
    in order to generate random pairs */
    shuffle(span4SDsinGBx.begin(), span4SDsinGBx.end(), gen);

    /* collide all randomly generated pairs of SDs */
    for (int i = 1; i < nsupers; i += 2)
    {
      collide_superdroplet_pair(gen, span4SDsinGBx[i-1].superdrop,
                                span4SDsinGBx[i].superdrop, scale_p, VOLUME);
    }
 
  }


  void collide_superdroplet_pair(std::mt19937 &gen, Superdrop &dropA,
                                 Superdrop &dropB, const int scale_p,
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
    const size_t eps1 = drop1.eps;
    const size_t eps2 = drop2.eps;

    /* 2. determine scaled probability of pair coalescence
    according to Shima et al. 2009 ("p_alpha" in paper) */
    const double prob_jk = pair_coalesce_probability(drop1, drop2, DELT, VOLUME);
    const double prob = scale_p * std::max(eps1, eps2) * prob_jk;

    /* 3. Monte Carlo step: randomly determine coalescence gamma factor */
    const size_t gamma = monte_carlo_gamma(gen, prob, eps1, eps2);

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

  size_t monte_carlo_gamma(std::mt19937 &gen, const double prob,
                        const size_t eps1, const size_t eps2) const
  /* calculates value of gamma factor in
  Monte Carlo collision-coalescence process
  according to Shima et al. 2009 */
  {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    const double phi = dis(gen); // random number phi in range [0,1]

    size_t gamma = 0;
    if (phi < (prob - floor(prob)))
    {
      gamma = floor(prob) + 1;
    }
    else if (phi >= (prob - floor(prob)))
    {
      gamma = floor(prob);
    }

    const size_t maxgamma = floor((eps1) / (eps2));

    std::cout << "gamma = " << gamma << " or " << maxgamma << "\n";
    std::cout << "floor = " <<  (eps1) << ", " << (eps2) << "\n";
  
    return std::min(gamma, maxgamma);
  }

  void coalesce_superdroplet_pair(Superdrop &drop1, Superdrop &drop2,
                                  const size_t gamma) const
  /* coalesce pair of superdroplets by changing multiplicity,
  radius and solute mass of each superdroplet in pair
  according to Shima et al. 2009 Section 5.1.3. part (5) */
  {
    if (drop1.eps == gamma * (drop2.eps))
    {
      twin_superdroplet_coalescence(drop1, drop2, gamma);
    }

    else if ((drop1.eps) > gamma * (drop2.eps))
    {
      different_superdroplet_coalescence(drop1, drop2, gamma);
    }

    else
    {
      std::string errormsg = "something undefined occured during colllision-coalescence" +
                             std::to_string(drop1.eps) + " < " +
                             //std::to_string(gamma * (drop2.eps)) + " ?";
                             std::to_string(gamma * (drop2.eps)) + ", " +
                             std::to_string(gamma) + ", " +
                             std::to_string(drop2.eps) + "?\n";
      throw std::invalid_argument(errormsg);
    }
  }

  void twin_superdroplet_coalescence(Superdrop &drop1,
                                     Superdrop &drop2,
                                     const size_t gamma) const
  /* if eps1 = gamma*eps2 coalescence makes twin SDs
  with same eps, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */
  { 
    const size_t new_eps = (drop2.eps) / 2.0;
    const double new_m_sol = (drop2.m_sol) + gamma * (drop1.m_sol);
    const double new_rcubed = pow((drop2.radius), 3.0) + gamma * (pow((drop1.radius), 3.0));
    const double new_r = pow(new_rcubed, 1.0 / 3.0);

    drop1.eps = new_eps;
    drop2.eps = (drop2.eps) - new_eps;

    drop1.radius = new_r;
    drop2.radius = new_r;

    drop1.m_sol = new_m_sol;
    drop2.m_sol = new_m_sol;
  }

  void different_superdroplet_coalescence(Superdrop &drop1,
                                          Superdrop &drop2,
                                          const size_t gamma) const
  /* if eps1 > gamma*eps2 coalescence grows drop2 radius and mass
  via decreasing multiplicity of drop1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */
  {
    drop1.eps = (drop1.eps) - gamma * (drop2.eps);

    const double new_rcubed = pow((drop2.radius), 3.0) + gamma * (pow((drop1.radius), 3.0));
    drop2.radius = pow(new_rcubed, 1.0 / 3.0);
    drop2.m_sol = (drop2.m_sol) + gamma * (drop1.m_sol);
  }

public:
  CollisionsMethod(const double DELT, PairCoalescenceProbability p)
      : DELT(DELT), pair_coalesce_probability(p) {}

  inline void operator()(const int currenttimestep,
                  std::span<SuperdropWithGridbox> span4SDsinGBx,
                  ThermoState &state,
                  std::mt19937 &gen) const
  /* this operator is used as an "adaptor" for using a run_step
  function in order to call collide_superdroplets. (*hint* run_step
  is usually found within a type that satisfies the SdmProcess concept) */
  {
    const double VOLUME = state.volume * pow(dlc::COORD0, 3.0); // volume in which collisions occur [m^3]
    collide_superdroplets(span4SDsinGBx, gen, VOLUME);
  }
};

#endif // COLLISIONSMETHOD_HPP