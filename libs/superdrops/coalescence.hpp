/*
 * ----- CLEO -----
 * File: coalescence.hpp
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
 * class and function to enact collision-coalescence events
 * in superdroplet model according to Shima et al. 2009.
 * Coalescence struct satisfies SDinGBxPairEnactX concept
 * used in Collisions struct */

#ifndef COALESCENCE_HPP
#define COALESCENCE_HPP

#include <cassert>

#include <Kokkos_Core.hpp>

#include "./collisions.hpp"
#include "./kokkosaliases_sd.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"

struct DoCoalescence
{
private:
  KOKKOS_FUNCTION unsigned long long
  coalescence_gamma(const unsigned long long xi1,
                    const unsigned long long xi2,
                    const double prob,
                    const double phi) const;
  /* calculates value of gamma factor in Monte Carlo
  collision-coalescence as in Shima et al. 2009 */

  KOKKOS_FUNCTION void
  coalesce_superdroplet_pair(Superdrop &drop1,
                             Superdrop &drop2,
                             const unsigned long long gamma) const;
  /* coalesce pair of superdroplets by changing multiplicity,
  radius and solute mass of each superdroplet in pair
  according to Shima et al. 2009 Section 5.1.3. part (5) */

  KOKKOS_FUNCTION void
  twin_superdroplet_coalescence(Superdrop &drop1,
                                Superdrop &drop2,
                                const unsigned long long gamma) const;
  /* if eps1 = gamma*eps2 coalescence makes twin SDs
  with same eps, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */

  KOKKOS_FUNCTION void
  different_superdroplet_coalescence(Superdrop &sd1,
                                     Superdrop &sd2,
                                     const unsigned long long gamma) const;
  /* if eps1 > gamma*eps2 coalescence grows sd2 radius and mass
  via decreasing multiplicity of sd1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */

public:
  KOKKOS_FUNCTION
  void operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const;
  /* this operator is used as an "adaptor" for using
  DoCoalescence as a function in DoCollisions that
  satistfies the PairEnactX concept */
};

template <PairProbability Probability>
inline MicrophysicalProcess auto
CollCoal(const unsigned int interval,
         const std::function<double(int)> int2realtime,
         const Probability collcoalprob)
/* constructs Microphysical Process for collision-coalescence
of superdroplets with a constant timestep 'interval' and
probability of collision-coalescence determined by 'collcoalprob' */
{
  const double DELT(int2realtime(interval));

  const DoCoalescence coal{};

  const DoCollisions<Probability, DoCoalescence> colls(DELT,
                                                       collcoalprob,
                                                       coal);
  return ConstTstepMicrophysics(interval, colls);
}

#endif // COLLISIONCOALESCENCE_HPP