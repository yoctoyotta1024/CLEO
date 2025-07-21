/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: breakup.hpp
 * Project: collisions
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality to enact collision-breakup events
 * in SDM analagous to to Shima et al. 2009.
 * Breakup struct satisfies PairEnactX
 * concept used in Collisions struct
 */

#ifndef LIBS_SUPERDROPS_COLLISIONS_BREAKUP_HPP_
#define LIBS_SUPERDROPS_COLLISIONS_BREAKUP_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <functional>

#include "../microphysicalprocess.hpp"
#include "../superdrop.hpp"
#include "./breakup_nfrags.hpp"
#include "./collisions.hpp"

template <NFragments NFrags>
struct DoBreakup {
 private:
  NFrags nfrags;

  /* calculates value of gamma factor in Monte Carlo
  collision-breakup, adapted from gamma for collision-
  coalescence in Shima et al. 2009. Here is is assumed
  maximally 1 breakup event can occur (gamma = 0 or 1)
  irrespective of whether scaled probability, prob > 1 */
  KOKKOS_FUNCTION
  unsigned int breakup_gamma(const double prob, const double phi) const;

  /* if xi1 = gamma*xi2 breakup of same multiplicity
  superdroplets produces (non-identical) twin superdroplets.
  Similar to Shima et al. 2009 Section 5.1.3. part (5) option (b).
  Note implicit assumption that gamma factor = 1. */
  KOKKOS_FUNCTION void twin_superdroplet_breakup(Superdrop &drop1, Superdrop &drop2) const;

  /* if xi1 > gamma*xi2 breakup alters drop2 radius and
  mass via decreasing multiplicity of drop1. Similar to
  Shima et al. 2009 Section 5.1.3. part (5) option (a).
  Note implicit assumption that gamma factor = 1. */
  KOKKOS_FUNCTION void different_superdroplet_breakup(Superdrop &drop1, Superdrop &drop2) const;

 public:
  explicit DoBreakup(const NFrags nfrags) : nfrags(nfrags) {}

  /**
   * @brief Operator used as an adaptor such that DoBreakup satisfies the PairEnactX concept
   * and so can be used as the EnactCollision function-like object in the DoCollisions struct.
   *
   * This operator calls functions to enact the collision-breakup of two super-droplets.
   *
   * @param drop1 First superdroplet.
   * @param drop2 Second superdroplet.
   * @param prob Probability of collision.
   * @param phi Phi value.
   * @return True if the resulting superdroplet is null, otherwise false.
   */
  KOKKOS_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2, const double prob, const double phi) const;

  /* enact collisional-breakup of droplets by changing
  multiplicity, radius and solute mass of each
  superdroplet in a pair. Method created by Author
  (no citation yet available). Note implicit assumption
  that gamma factor = 1. */
  KOKKOS_FUNCTION void breakup_superdroplet_pair(Superdrop &drop1, Superdrop &drop2) const;
};

/* constructs Microphysical Process for collision-breakup
of superdroplets with a constant timestep 'interval' and
probability of collision-breakup determined by 'collbuprob' */
template <PairProbability Probability, NFragments NFrags>
inline MicrophysicalProcess auto CollBu(const unsigned int interval,
                                        const std::function<double(unsigned int)> int2realtime,
                                        const Probability collbuprob, const NFrags nfrags) {
  const auto DELT = double{int2realtime(interval)};

  const DoBreakup bu(nfrags);
  const MicrophysicsFunc auto colls =
      DoCollisions<Probability, DoBreakup<NFrags>>(DELT, collbuprob, bu);

  return ConstTstepMicrophysics(interval, colls);
}

/**
 * @brief Operator used as an adaptor such that DoBreakup satisfies the PairEnactX concept
 * and so can be used as the EnactCollision function-like object in the DoCollisions struct.
 *
 * This operator calls functions to enact the collision-breakup of two super-droplets.
 *
 * @param drop1 First superdroplet.
 * @param drop2 Second superdroplet.
 * @param prob Probability of collision.
 * @param phi Phi value.
 * @return True if the resulting superdroplet is null, otherwise false.
 */
template <NFragments NFrags>
KOKKOS_FUNCTION bool DoBreakup<NFrags>::operator()(Superdrop &drop1, Superdrop &drop2,
                                                   const double prob, const double phi) const {
  /* enact collision-breakup on pair of superdroplets if
  gamma factor for collision-breakup is not zero */
  if (breakup_gamma(prob, phi) != 0) {
    breakup_superdroplet_pair(drop1, drop2);
  }

  return 0;
}

/* calculates value of gamma factor in Monte Carlo
collision-breakup, adapted from gamma for collision-
coalescence in Shima et al. 2009. Here is is assumed
maximally 1 breakup event can occur (gamma = 0 or 1)
irrespective of whether scaled probability, prob > 1 */
template <NFragments NFrags>
KOKKOS_FUNCTION unsigned int DoBreakup<NFrags>::breakup_gamma(const double prob,
                                                              const double phi) const {
  if (phi < (prob - Kokkos::floor(prob))) {
    return 1;
  } else {  // if phi >= (prob - floor(prob))
    return 0;
  }
}

/* enact collisional-breakup of droplets by changing
multiplicity, radius and solute mass of each
superdroplet in a pair. Method created by Author
(no citation yet available). Note implicit assumption
that gamma factor = 1. */
template <NFragments NFrags>
KOKKOS_FUNCTION void DoBreakup<NFrags>::breakup_superdroplet_pair(Superdrop &drop1,
                                                                  Superdrop &drop2) const {
  if (drop1.get_xi() == drop2.get_xi()) {
    twin_superdroplet_breakup(drop1, drop2);
  } else {
    different_superdroplet_breakup(drop1, drop2);
  }
}

/* if xi1 = gamma*xi2 breakup of same multiplicity
superdroplets produces (non-identical) twin superdroplets.
Similar to Shima et al. 2009 Section 5.1.3. part (5) option (b).
Note implicit assumption that gamma factor = 1.
_Note:_ Implicit casting of xi from uint64_t to double. */
template <NFragments NFrags>
KOKKOS_FUNCTION void DoBreakup<NFrags>::twin_superdroplet_breakup(Superdrop &drop1,
                                                                  Superdrop &drop2) const {
  const auto old_xi = drop2.get_xi();  // = drop1.xi
  const auto totnfrags = double{nfrags(drop1, drop2) * old_xi};
  assert(((totnfrags / old_xi) > 2.5) && "nfrags must be > 2.5");

  const auto new_xi1 = static_cast<uint64_t>(Kokkos::round(totnfrags / 2));
  const auto new_xi2 = static_cast<uint64_t>(Kokkos::round(totnfrags - new_xi1));
  const auto new_xitot = new_xi1 + new_xi2;
  assert((new_xi2 > old_xi) && "nfrags must increase the drop2's multiplicity during breakup");
  assert((new_xitot > (old_xi * 2)) && "nfrags must increase total multiplicity during breakup");

  const auto sum_rcubed = double{drop1.rcubed() + drop2.rcubed()};
  const auto new_rcubed = double{sum_rcubed * old_xi / new_xitot};
  const auto new_r = double{Kokkos::pow(new_rcubed, (1.0 / 3.0))};

  const auto sum_msol = double{drop1.get_msol() + drop2.get_msol()};
  const auto new_msol = double{sum_msol * old_xi / new_xitot};

  drop1.set_xi(new_xi1);
  drop2.set_xi(new_xi2);

  drop1.set_radius(new_r);
  drop2.set_radius(new_r);

  drop1.set_msol(new_msol);
  drop2.set_msol(new_msol);
}

/* if xi1 > gamma*xi2 breakup alters drop2 radius and
mass via decreasing multiplicity of drop1. Similar to
Shima et al. 2009 Section 5.1.3. part (5) option (a).
Note implicit assumption that gamma factor = 1.
_Note:_ Implicit casting of xi from uint64_t to double. */
template <NFragments NFrags>
KOKKOS_FUNCTION void DoBreakup<NFrags>::different_superdroplet_breakup(Superdrop &drop1,
                                                                       Superdrop &drop2) const {
  const auto old_xi1 = drop1.get_xi();
  const auto old_xi2 = drop2.get_xi();

  const auto new_xi1 = old_xi1 - old_xi2;
  drop1.set_xi(new_xi1);

  const auto totnfrags = double{nfrags(drop1, drop2) * old_xi2};
  const auto new_xi2 = static_cast<uint64_t>(Kokkos::round(totnfrags));
  assert(((totnfrags / old_xi2) > 2.5) && "nfrags must be > 2.5");

  assert((new_xi2 > old_xi2) && "nfrags must increase the drop2's multiplicity during breakup");
  assert(((new_xi1 + new_xi2) > (old_xi1 + old_xi2)) &&
         "nfrags must increase the total multiplicity during breakup");

  const auto sum_rcubed = double{drop1.rcubed() + drop2.rcubed()};
  const auto new_rcubed = double{sum_rcubed * old_xi2 / new_xi2};
  const auto new_r = double{Kokkos::pow(new_rcubed, (1.0 / 3.0))};

  const auto sum_msol = double{drop1.get_msol() + drop2.get_msol()};
  const auto new_msol = double{sum_msol * old_xi2 / new_xi2};

  drop2.set_xi(new_xi2);
  drop2.set_radius(new_r);
  drop2.set_msol(new_msol);
}

#endif  // LIBS_SUPERDROPS_COLLISIONS_BREAKUP_HPP_
