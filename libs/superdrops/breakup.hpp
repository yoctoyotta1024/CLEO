/*
 * ----- CLEO -----
 * File: breakup.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 14th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality to enact collision-breakup events
 * in SDM analagous to to Shima et al. 2009.
 * Breakup struct satisfies PairEnactX
 * concept used in Collisions struct
 */

#ifndef BREAKUP_HPP
#define BREAKUP_HPP

#include <functional>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "./collisions.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"
#include "./breakup_nfrags.hpp"

template <NFragments NFrags>
struct DoBreakup
{
private:
  NFrags nfrags;

  KOKKOS_FUNCTION
  unsigned int breakup_gamma(const double prob,
                             const double phi) const;
  /* calculates value of gamma factor in Monte Carlo
  collision-breakup, adapted from gamma for collision-
  coalescence in Shima et al. 2009. Here is is assumed
  maximally 1 breakup event can occur (gamma = 0 or 1)
  irrespective of whether scaled probability, prob > 1 */

  KOKKOS_FUNCTION void
  twin_superdroplet_breakup(Superdrop &drop1,
                            Superdrop &drop2) const;
  /* if xi1 = gamma*xi2 breakup of same multiplicity
  superdroplets produces (non-identical) twin superdroplets.
  Similar to Shima et al. 2009 Section 5.1.3. part (5) option (b).
  Note implicit assumption that gamma factor = 1. */

  KOKKOS_FUNCTION void
  different_superdroplet_breakup(Superdrop &drop1,
                                 Superdrop &drop2) const;
  /* if xi1 > gamma*xi2 breakup alters drop2 radius and
  mass via decreasing multiplicity of drop1. Similar to
  Shima et al. 2009 Section 5.1.3. part (5) option (a).
  Note implicit assumption that gamma factor = 1. */

public:
  DoBreakup(const NFrags nfrags) : nfrags(nfrags) {}

  KOKKOS_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const;
  /* this operator is used as an "adaptor" for
  using DoBreakup as a function in DoCollisions
  that satistfies the PairEnactX concept */

  KOKKOS_FUNCTION void
  breakup_superdroplet_pair(Superdrop &drop1,
                            Superdrop &drop2) const;
  /* enact collisional-breakup of droplets by changing
  multiplicity, radius and solute mass of each
  superdroplet in a pair. Method created by Author
  (no citation yet available). Note implicit assumption
  that gamma factor = 1. */
};

template <PairProbability Probability, NFragments NFrags>
inline MicrophysicalProcess auto
CollBu(const unsigned int interval,
       const std::function<double(unsigned int)> int2realtime,
       const Probability collbuprob,
       const NFrags nfrags)
/* constructs Microphysical Process for collision-breakup
of superdroplets with a constant timestep 'interval' and
probability of collision-breakup determined by 'collbuprob' */
{
  const auto DELT = double{int2realtime(interval)};

  const DoBreakup bu(nfrags);
  const DoCollisions<Probability, DoBreakup<NFrags>> colls(DELT,
                                                           collbuprob,
                                                           bu);

  return ConstTstepMicrophysics(interval, colls);
}

template <NFragments NFrags>
KOKKOS_FUNCTION bool
DoBreakup<NFrags>::operator()(Superdrop &drop1, Superdrop &drop2,
                              const double prob, const double phi) const
/* this operator is used as an "adaptor" for
using DoBreakup as a function in DoCollisions
that satistfies the PairEnactX concept */
{
  /* enact collision-breakup on pair of superdroplets if
  gamma factor for collision-breakup is not zero */
  if (breakup_gamma(prob, phi) != 0)
  {
    breakup_superdroplet_pair(drop1, drop2);
  }

  return 0;
}

template <NFragments NFrags>
KOKKOS_FUNCTION unsigned int
DoBreakup<NFrags>::breakup_gamma(const double prob,
                                 const double phi) const
/* calculates value of gamma factor in Monte Carlo
collision-breakup, adapted from gamma for collision-
coalescence in Shima et al. 2009. Here is is assumed
maximally 1 breakup event can occur (gamma = 0 or 1)
irrespective of whether scaled probability, prob > 1 */
{
  if (phi < (prob - floor(prob)))
  {
    return 1;
  }
  else // if phi >= (prob - floor(prob))
  {
    return 0;
  }
}

template <NFragments NFrags>
KOKKOS_FUNCTION void
DoBreakup<NFrags>::breakup_superdroplet_pair(Superdrop &drop1,
                                             Superdrop &drop2) const
/* enact collisional-breakup of droplets by changing
multiplicity, radius and solute mass of each
superdroplet in a pair. Method created by Author
(no citation yet available). Note implicit assumption
that gamma factor = 1. */
{
  if (drop1.get_xi() == drop2.get_xi())
  {
    twin_superdroplet_breakup(drop1, drop2);
  }

  else
  {
    different_superdroplet_breakup(drop1, drop2);
  }
}

template <NFragments NFrags>
KOKKOS_FUNCTION void
DoBreakup<NFrags>::twin_superdroplet_breakup(Superdrop &drop1,
                                             Superdrop &drop2) const
/* if xi1 = gamma*xi2 breakup of same multiplicity
superdroplets produces (non-identical) twin superdroplets.
Similar to Shima et al. 2009 Section 5.1.3. part (5) option (b).
Note implicit assumption that gamma factor = 1.
Note: implicit casting of xi from unsigned long long to double. */
{
  const auto old_xi = drop2.get_xi(); // = drop1.xi
  const auto totnfrags = double{nfrags(drop1, drop2) * old_xi};
  const auto new_xi = (unsigned long long)Kokkos::round(totnfrags) / 2;

  const auto sum_rcubed = double{drop1.rcubed() + drop2.rcubed()};
  const auto new_rcubed = double{sum_rcubed * old_xi / new_xi};
  const auto new_r = double{Kokkos::pow(new_rcubed, (1.0 / 3.0))};

  const auto sum_msol = double{drop1.get_msol() + drop2.get_msol()};
  const auto new_msol = double{sum_msol * old_xi / new_xi};

  drop1.set_xi(new_xi);
  drop2.set_xi(old_xi - new_xi);

  drop1.set_radius(new_r);
  drop2.set_radius(new_r);

  drop1.set_msol(new_msol);
  drop2.set_msol(new_msol);
}

template <NFragments NFrags>
KOKKOS_FUNCTION void
DoBreakup<NFrags>::different_superdroplet_breakup(Superdrop &drop1,
                                                  Superdrop &drop2) const
/* if xi1 > gamma*xi2 breakup alters drop2 radius and
mass via decreasing multiplicity of drop1. Similar to
Shima et al. 2009 Section 5.1.3. part (5) option (a).
Note implicit assumption that gamma factor = 1.
Note: implicit casting of xi from unsigned long long to double. */
{
  drop1.set_xi(drop1.get_xi() - drop2.get_xi());

  const auto old_xi = drop2.get_xi();
  const auto totnfrags = double{nfrags(drop1, drop2) * old_xi};
  const auto new_xi = (unsigned long long)Kokkos::round(totnfrags);

  const auto sum_rcubed = double{drop1.rcubed() + drop2.rcubed()};
  const auto new_rcubed = double{sum_rcubed * old_xi / new_xi};
  const auto new_r = double{Kokkos::pow(new_rcubed, (1.0 / 3.0))};

  const auto sum_msol = double{drop1.get_msol() + drop2.get_msol()};
  const auto new_msol = double{sum_msol * old_xi / new_xi};

  drop2.set_xi(new_xi);
  drop2.set_radius(new_r); 
  drop2.set_msol(new_msol);
}

#endif // BREAKUP_HPP