/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: coalescence.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 11th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * class and function to enact collision-coalescence events
 * in superdroplet model according to Shima et al. 2009.
 * Coalescence struct satisfies PairEnactX concept
 * used in Collisions struct
 */

#ifndef LIBS_SUPERDROPS_COALESCENCE_HPP_
#define LIBS_SUPERDROPS_COALESCENCE_HPP_

#include <cassert>
#include <functional>

#include <Kokkos_Core.hpp>

#include "./collisions.hpp"
#include "./microphysicalprocess.hpp"
#include "./nullsuperdrops.hpp"
#include "./superdrop.hpp"

struct DoCoalescence {
 private:
  /* if xi1 = gamma*xi2 coalescence makes twin SDs
  with same xi, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */
  KOKKOS_FUNCTION void twin_superdroplet_coalescence(const uint64_t gamma, Superdrop &drop1,
                                                     Superdrop &drop2) const;

  /* if xi1 > gamma*xi2 coalescence grows sd2 radius and mass
  via decreasing multiplicity of sd1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */
  KOKKOS_FUNCTION void different_superdroplet_coalescence(const uint64_t gamma,
                                                          Superdrop &drop1,
                                                          Superdrop &drop2) const;

 public:
  /* this operator is used as an "adaptor" for using
  DoCoalescence as a function in DoCollisions that
  satistfies the PairEnactX concept */
  KOKKOS_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2, const double prob, const double phi) const;

  /* calculates value of gamma factor in Monte Carlo
  collision-coalescence as in Shima et al. 2009 */
  KOKKOS_FUNCTION uint64_t coalescence_gamma(const uint64_t xi1, const uint64_t xi2,
                                             const double prob, const double phi) const;

  /* coalesce pair of superdroplets by changing multiplicity,
  radius and solute mass of each superdroplet in pair
  according to Shima et al. 2009 Section 5.1.3. part (5) */
  KOKKOS_FUNCTION bool coalesce_superdroplet_pair(const uint64_t gamma, Superdrop &drop1,
                                                  Superdrop &drop2) const;
};

/* constructs Microphysical Process for collision-coalescence
of superdroplets with a constant timestep 'interval' and
probability of collision-coalescence determined by 'collcoalprob' */
template <PairProbability Probability>
inline MicrophysicalProcess auto CollCoal(const unsigned int interval,
                                          const std::function<double(unsigned int)> int2realtime,
                                          const Probability collcoalprob) {
  const auto DELT = int2realtime(interval);

  const DoCoalescence coal{};
  const DoCollisions<Probability, DoCoalescence> colls(DELT, collcoalprob, coal);
  return ConstTstepMicrophysics(interval, colls);
}

#endif   // LIBS_SUPERDROPS_COALESCENCE_HPP_
