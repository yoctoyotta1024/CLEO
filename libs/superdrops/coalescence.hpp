/*
 * ----- CLEO -----
 * File: coalescence.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 22nd December 2023
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
 * Coalescence struct satisfies PairEnactX concept
 * used in Collisions struct */

#ifndef COALESCENCE_HPP
#define COALESCENCE_HPP

#include <cassert>
#include <functional>

#include <Kokkos_Core.hpp>

#include "./collisions.hpp"
#include "./microphysicalprocess.hpp"
#include "./nullsuperdrops.hpp"
#include "./superdrop.hpp"

struct DoCoalescence
{
private:
  KOKKOS_FUNCTION void
  twin_superdroplet_coalescence(const unsigned long long gamma,
                                Superdrop &drop1,
                                Superdrop &drop2) const;
  /* if xi1 = gamma*xi2 coalescence makes twin SDs
  with same xi, r and solute mass. According to Shima et al. 2009
  Section 5.1.3. part (5) option (b)  */

  KOKKOS_FUNCTION void
  different_superdroplet_coalescence(const unsigned long long gamma,
                                     Superdrop &drop1,
                                     Superdrop &drop2) const;
  /* if xi1 > gamma*xi2 coalescence grows sd2 radius and mass
  via decreasing multiplicity of sd1. According to
  Shima et al. 2009 Section 5.1.3. part (5) option (a)  */

public:
  KOKKOS_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const;
  /* this operator is used as an "adaptor" for using
  DoCoalescence as a function in DoCollisions that
  satistfies the PairEnactX concept */

  KOKKOS_FUNCTION unsigned long long
  coalescence_gamma(const unsigned long long xi1,
                    const unsigned long long xi2,
                    const double prob,
                    const double phi) const;
  /* calculates value of gamma factor in Monte Carlo
  collision-coalescence as in Shima et al. 2009 */

  KOKKOS_FUNCTION bool
  coalesce_superdroplet_pair(const unsigned long long gamma,
                             Superdrop &drop1,
                             Superdrop &drop2) const;
  /* coalesce pair of superdroplets by changing multiplicity,
  radius and solute mass of each superdroplet in pair
  according to Shima et al. 2009 Section 5.1.3. part (5) */
};

template <PairProbability Probability>
inline MicrophysicalProcess auto
CollCoal(const unsigned int interval,
         const std::function<double(unsigned int)> int2realtime,
         const Probability collcoalprob)
/* constructs Microphysical Process for collision-coalescence
of superdroplets with a constant timestep 'interval' and
probability of collision-coalescence determined by 'collcoalprob' */
{
  const auto DELT = int2realtime(interval);

  const DoCoalescence coal{};
  const DoCollisions<Probability, DoCoalescence> colls(DELT,
                                                       collcoalprob,
                                                       coal);
  return ConstTstepMicrophysics(interval, colls);
}

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

KOKKOS_FUNCTION bool
DoCoalescence::operator()(Superdrop &drop1, Superdrop &drop2,
                          const double prob, const double phi) const
/* this operator is used as an "adaptor" for using
DoCoalescence as a function in DoCollisions that
satistfies the PairEnactX concept */
{
  /* 1. calculate gamma factor for collision-coalescence  */
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();
  const auto gamma = coalescence_gamma(xi1, xi2, prob, phi);

  /* 2. enact collision-coalescence on pair
  of superdroplets if gamma is not zero */
  if (gamma != 0)
  {
    return coalesce_superdroplet_pair(gamma, drop1, drop2);
  }

  return 0;
}

KOKKOS_FUNCTION unsigned long long
DoCoalescence::coalescence_gamma(const unsigned long long xi1,
                                 const unsigned long long xi2,
                                 const double prob,
                                 const double phi) const
/* calculates value of gamma factor in Monte Carlo
collision-coalescence as in Shima et al. 2009 */
{
  unsigned long long gamma = floor(prob); // if phi >= (prob - floor(prob))
  if (phi < (prob - gamma))
  {
    ++gamma;
  }

  const auto maxgamma = xi1 / xi2; // same as floor() for positive ints

  return Kokkos::fmin(gamma, maxgamma);
}

KOKKOS_FUNCTION bool
DoCoalescence::coalesce_superdroplet_pair(const unsigned long long gamma,
                                          Superdrop &drop1,
                                          Superdrop &drop2) const
/* coalesce pair of superdroplets by changing multiplicity,
radius and solute mass of each superdroplet in pair
according to Shima et al. 2009 Section 5.1.3. part (5) */
{
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();

  if (xi1 - gamma * xi2 > 0)
  {
    different_superdroplet_coalescence(gamma, drop1, drop2);
    return 0;
  }

  else if (xi1 - gamma * xi2 == 0)
  {
    twin_superdroplet_coalescence(gamma, drop1, drop2);
    
    /* if xi1 = xi2 = 1 before coalesence, then xi1=0 now */
    return is_null_superdrop(drop1);
    // return if_null_superdrop(drop1);
  }

  assert((xi1 >= gamma * xi2) && "something undefined occured "
                                   "during colllision-coalescence");
  return 0;                               
}

KOKKOS_FUNCTION void
DoCoalescence::twin_superdroplet_coalescence(const unsigned long long gamma,
                                             Superdrop &drop1,
                                             Superdrop &drop2) const
/* if xi1 = gamma*xi2 coalescence makes twin SDs
with same xi, r and solute mass. According to Shima et al. 2009
Section 5.1.3. part (5) option (b).
Note: implicit casting of gamma (i.e. therefore droplets'
xi values) from unsigned long long to double. */
{
  const auto old_xi = drop2.get_xi(); // = drop1.xi
  const auto new_xi = old_xi / 2; // same as floor() for positive ints

  const auto new_rcubed = double{drop2.rcubed() +
                                 gamma * drop1.rcubed()};
  const auto new_r = double{Kokkos::pow(new_rcubed, (1.0 / 3.0))};

  const auto new_msol = double{drop2.get_msol() +
                               gamma * drop1.get_msol()};

  drop1.set_xi(new_xi);
  drop2.set_xi(old_xi - new_xi);

  drop1.set_radius(new_r);
  drop2.set_radius(new_r);

  drop1.set_msol(new_msol);
  drop2.set_msol(new_msol);
}

KOKKOS_FUNCTION void
DoCoalescence::different_superdroplet_coalescence(const unsigned long long gamma,
                                                  Superdrop &drop1,
                                                  Superdrop &drop2) const
/* if xi1 > gamma*xi2 coalescence grows drop2 radius and mass
via decreasing multiplicity of drop1. According to
Shima et al. 2009 Section 5.1.3. part (5) option (a)
Note: implicit casting of gamma (i.e. therefore droplets'
xi values) from unsigned long long to double. */
{
  drop1.set_xi(drop1.get_xi() - gamma * drop2.get_xi());

  const auto new_rcubed = double{drop2.rcubed() +
                                 gamma * drop1.rcubed()};

  drop2.set_radius(Kokkos::pow(new_rcubed, (1.0 / 3.0)));
  drop2.set_msol(drop2.get_msol() + gamma * drop1.get_msol());
}

#endif // COLLISIONCOALESCENCE_HPP