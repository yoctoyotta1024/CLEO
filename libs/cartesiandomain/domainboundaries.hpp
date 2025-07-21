/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: domainboundaries.hpp
 * Project: cartesiandomain
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functions to implement finite or periodic boundary conditions
 * at the edges of a cartesian domain e.g. for returning the neighbouring
 * gbxindex and value of a superdroplet's coord when the superdroplet
 * crosses the domain boundary in a particular direction.
 */

#ifndef LIBS_CARTESIANDOMAIN_DOMAINBOUNDARIES_HPP_
#define LIBS_CARTESIANDOMAIN_DOMAINBOUNDARIES_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "../cleoconstants.hpp"

/* returns true if gbxindex for gridbox is at/beyond a
boundary of cartesian domain, given neighbouring indexes
are +/- increment from gbxindex and the number of
gridboxes making up the domain in that direction (ndim) */
KOKKOS_INLINE_FUNCTION bool beyond_domainboundary(const unsigned int gbxindex,
                                                  const unsigned int increment, const size_t ndim) {
  return (gbxindex / increment) % ndim == 0;
}

/* returns {backwards, forwards} gridbox neighbours with
treatment of neighbours as if bounds of domain are finite.
This means that no neighbour exists above highest / below lowest
gridboxes in a given direction. Non-existent neighbours are
defined with gbxindex = max unsigned int, meaning in a given
direction the gbxindex of the backwards / forwards neighbour
of a gridbox at the edge of the domain is a max unsigned int */
inline Kokkos::pair<unsigned int, unsigned int> finitedomain_nghbrs(const unsigned int idx,
                                                                    const unsigned int increment,
                                                                    const unsigned int ndim) {
  unsigned int forward(idx + increment);
  unsigned int backward(idx - increment);

  // at lower edge of domain
  if (beyond_domainboundary(idx, increment, ndim)) {
    backward = LIMITVALUES::oob_gbxindex;
  }

  // at upper edge of domain
  if (beyond_domainboundary(forward, increment, ndim)) {
    forward = LIMITVALUES::oob_gbxindex;
  }

  return {backward, forward};
}

/* returns {backwards, forwards} gridbox neighbours with
treatment of neighbours as if bounds of domain are periodic.
This means that highest and lowest gridboxes in a given
direction are each others' neighbours. ie. index of neighbour
forwards of gridboxes at the uppermost edge of domain is the
lowermost gridbox in that direction (and vice versa). */
inline Kokkos::pair<unsigned int, unsigned int> periodicdomain_nghbrs(const unsigned int idx,
                                                                      const unsigned int increment,
                                                                      const unsigned int ndim) {
  unsigned int forward(idx + increment);
  unsigned int backward(idx - increment);

  // at lower edge of domain
  if (beyond_domainboundary(idx, increment, ndim)) {
    backward = idx + (ndim - 1) * increment;
  }

  // at upper edge of domain
  if (beyond_domainboundary(forward, increment, ndim)) {
    forward = idx - (ndim - 1) * increment;
  }

  return {backward, forward};
}

/* Finite domain boundaries don't
change superdroplet coord */
KOKKOS_INLINE_FUNCTION double coordbeyond_finitedomain(const double coord, const double lim1,
                                                       const double lim2) {
  return coord;
}

/* In periodic domain, two scenarios:
a) If superdroplet crosses lower boundary of domain,
lim1 = upper bound of backwards neighbour from gbx (upper boundary of domain)
lim2 = lower bound of gridbox (lower boundary of domain) so
coord -> coord + length_of_domain
b) If superdroplet crosses upper boundary of domain,
lim1 = lower bound of forwards neighbour from gbx (lower boundary of domain)
lim2 = upper bound of gridbox (upper boundary of domain) so
coord -> coord - length_of_domain */
KOKKOS_INLINE_FUNCTION double coordbeyond_periodicdomain(const double coord, const double lim1,
                                                         const double lim2) {
  return coord + lim1 - lim2;  // periodic domain coord -> coord +/- |length_of_domain|
}

#endif  // LIBS_CARTESIANDOMAIN_DOMAINBOUNDARIES_HPP_
