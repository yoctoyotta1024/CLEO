// Author: Clara Bayley
// File: coal_breakup.hpp
/* Header file for class that enacts
collision events in which either coalescence,
or breakup occurs with fixed efficiency
Ecoal (and 1 - Ecoal). ConstCoalBreakup
struct satisfies SDPairEnactX concept used in
CollisionX struct */

#ifndef COAL_BREAKUP_HPP
#define COAL_BREAKUP_HPP

#include <functional>
#include <concepts>
#include <stdexcept>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"
#include "./collisionxkernels.hpp"
#include "./collisionx.hpp"
#include "./coalescence.hpp"
#include "./breakup.hpp"

class CoalBreakupConstEff
/* class is method for coalescence / breakup between
two superdroplets with constant efficiency of
coalescence coaleff. (Can be used in collisionsx
struct to enact collision-coalescence or
collision-breakup events in the superdroplet model) */
{
private:
  Coalescence coal;
  Breakup breakup;
  double coaleff;
  double bueff;

public:
  CoalBreakupConstEff(const double infrags, const double coaleff)
      : coal(Coalescence{}), breakup(Breakup(infrags)),
        coaleff(coaleff), bueff(1.0 - coaleff)
  {
    if ((coaleff > 1.0) || (coaleff < 0.0))
    {
      throw std::invalid_argument("Invalid coalescence efficiency, coaleff");
    }
  }
  
  void operator()(SuperdropWithGbxindex &SDinGBx1,
                  SuperdropWithGbxindex &SDinGBx2,
                  const double probcoll, const double phi) const
  /* this operator is used as an "adaptor" for using
  CoalBreakupConstEff as a function in CollisionsX
  that satistfies the SDPairEnactX concept.
  *note* operator uses probcoll, probability of collision,
  NOT probability of collision-coalescence! */
  {
    const unsigned long long eps1(SDinGBx1.superdrop.eps);
    const unsigned long long eps2(SDinGBx2.superdrop.eps);

    /* 1. calculate gamma factor for collision-coalescence  */
    const double probcoal(probcoll * coaleff);
    const unsigned long long gamma_coal(coal.coalescence_gamma(eps1,
                                                               eps2,
                                                               probcoal,
                                                               phi));
    /* 2. enact collision-coalescence between pair
      of superdroplets if gamma is not zero */
    if (gamma_coal != 0)
    {
      coal.coalesce_superdroplet_pair(SDinGBx1, SDinGBx2, gamma_coal);
    }

    else // if not coalescence, check for breakup
    {
      /* 3. calculate gamma factor for collision-breakup  */
      const double probbu(probcoll * bueff);
      const unsigned long long gamma_bu(breakup.breakup_gamma(eps1,
                                                              eps2,
                                                              probbu,
                                                              phi));
      /* 4. enact collision-breakup between pair
        of superdroplets if gamma is not zero */
      if (gamma_bu != 0)
      {
        breakup.breakup_superdroplet_pair(SDinGBx1.superdrop,
                                          SDinGBx1.superdrop);
      }
    }
  }
};

SdmProcess auto
CollisionCoalBuConst(const int interval,
                     const std::function<double(int)> int2time,
                     const double nfrags,
                     const double coalrate,
                     const double burate)
/* SDM process for collisions of superdroplets
followed by coalescence or breakup with constant
coaleff (similar to de Jong et al. 2023 sect. 3) */
{
  const double realtstep = int2time(interval);

  const double kernel(coalrate + burate);
  const CollConstProb collprob(kernel);

  const double coaleff(coalrate / (coalrate + burate));
  CollisionX<CollConstProb, CoalBreakupConstEff>
      coalbu(realtstep, collprob,
             CoalBreakupConstEff(nfrags, coaleff));

  return ConstTstepProcess{interval, coalbu};
};

#endif // COAL_BREAKUP_HPP