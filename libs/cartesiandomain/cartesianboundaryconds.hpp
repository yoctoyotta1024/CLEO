/*
 * ----- CLEO -----
 * File: cartesianboundaryconds.hpp
 * Project: cartesiandomain
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 19th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions that determine the boundary conditions
 * at the edges of the cartesian domain e.g. for
 * returning the neighbouring gbxindex and
 * value of a superdroplet's coord when the superdroplet
 * crosses the domain boundary in a particular direction
 */

#ifndef CARTESIANBOUNDARYCONDS_HPP
#define CARTESIANBOUNDARYCONDS_HPP

/* NOTE: boundary conditions of domain are defined as:
  z: FINITE   (see cartesian_coord3nghbrs & coord3_beyondz)
  x: PERIODIC (see cartesian_coord3nghbrs & coord1_beyondx)
  y: PERIODIC (see cartesian_coord3nghbrs & coord2_beyondy)
*/

#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "../cleoconstants.hpp"

KOKKOS_INLINE_FUNCTION unsigned int
outofbounds_gbxindex()
/* value to set sdgbxindex to indicate superdrop
is out of domain (ie. not a valid gbxindex) */
{
  return LIMITVALUES::uintmax;
}

KOKKOS_INLINE_FUNCTION bool
beyond_domainboundary(const unsigned int gbxindex,
                           const unsigned int increment,
                           const size_t ndim)
/* returns true if gbxindex for gridbox is at/beyond a
boundary of cartesian domain, given neighbouring indexes
are +/- increment from gbxindex and the number of
gridboxes making up the domain in that direction (ndim) */
{
  return (gbxindex / increment) % ndim == 0;
}

inline Kokkos::pair<unsigned int, unsigned int>
finitedomain_nghbrs(const unsigned int idx,
                    const unsigned int increment,
                    const unsigned int ndim)
/* returns {backwards, forwards} gridbox neighbours with
treatment of neighbours as if bounds of domain are finite.
This means that no neighbour exists above highest / below lowest
gridboxes in a given direction. Non-existent neighbours are
defined with gbxindex = max unsigned int, meaning in a given
direction the gbxindex of the backwards / forwards neighbour
of a gridbox at the edge of the domain is a max unsigned int */
{
  unsigned int forward(idx + increment);
  unsigned int backward(idx - increment);

  if (beyond_domainboundary(idx, increment, ndim)) // at lower edge of domain
  {
    backward = outofbounds_gbxindex();
  }

  if (beyond_domainboundary(forward, increment, ndim)) // at upper edge of domain
  {
    forward = outofbounds_gbxindex();
  }

  return {backward, forward};
}

inline Kokkos::pair<unsigned int, unsigned int>
periodicdomain_nghbrs(const unsigned int idx,
                      const unsigned int increment,
                      const unsigned int ndim)
/* returns {backwards, forwards} gridbox neighbours with
treatment of neighbours as if bounds of domain are periodic.
This means that highest and lowest gridboxes in a given
direction are each others' neighbours. ie. index of neighbour
forwards of gridboxes at the uppermost edge of domain is the
lowermost gridbox in that direction (and vice versa). */
{
  unsigned int forward(idx + increment);
  unsigned int backward(idx - increment);

  if (beyond_domainboundary(idx, increment, ndim)) // at lower edge of domain
  {
    backward = idx + (ndim - 1) * increment;
  }

  if (beyond_domainboundary(forward, increment, ndim)) // at upper edge of domain
  {
    forward = idx - (ndim - 1) * increment;
  }

  return {backward, forward};
}

KOKKOS_INLINE_FUNCTION double
coordbeyond_finitedomain(const double coord,
                         const double lim1,
                         const double lim2)
/* Finite domain boundaries don't change superdroplet coord */
{
  return coord; // finite domain therefore don't change coord
}

KOKKOS_INLINE_FUNCTION double
coordbeyond_periodicdomain(const double coord,
                           const double lim1,
                           const double lim2)
/* In periodic domain, two scenarios:
a) If superdroplet crosses lower boundary of domain,
lim1 = upper bound of backwards neighbour from gbx (upper boundary of domain)
lim2 = lower bound of gridbox (lower boundary of domain) so
coord -> coord + length_of_domain
b) If superdroplet crosses upper boundary of domain,
lim1 = lower bound of forwards neighbour from gbx (lower boundary of domain)
lim2 = upper bound of gridbox (upper boundary of domain) so
coord -> coord - length_of_domain */
{
  return coord + lim1 - lim2; // periodic domain coord -> coord +/- |length_of_domain|
}

inline Kokkos::pair<unsigned int, unsigned int>
cartesian_coord3nghbrs(const unsigned int idx,
                  const std::vector<size_t> &ndims)
/* returns pair for gbx index of neighbour in the
{backwards, forwards} z direction given a gridbox with
gbxidx='idx' in a cartesian domain. Treatment of neighbours
for gridboxes at the edges of the domain is either finite
(null neighbour) or periodic (cyclic neighbour) */
{
  return finitedomain_nghbrs(idx, 1, ndims.at(0));
  // return periodicdomain_nghbrs(idx, 1, ndims.at(0));
}

inline Kokkos::pair<unsigned int, unsigned int>
cartesian_coord1nghbrs(const unsigned int idx,
                  const std::vector<size_t> &ndims)
/* returns pair for gbx index of neighbour in the
{backwards, forwards} x direction given a gridbox with
gbxidx='idx' in a cartesian domain. Treatment of neighbours
for gridboxes at the edges of the domain is either finite
(null neighbour) or periodic (cyclic neighbour) */
{
  const auto nz = (unsigned int)ndims.at(0); // no. gridboxes in z direction
  // return finitedomain_nghbrs(idx, nz, ndims.at(1));
  return periodicdomain_nghbrs(idx, nz, ndims.at(1));
}

inline Kokkos::pair<unsigned int, unsigned int>
cartesian_coord2nghbrs(const unsigned int idx,
                  const std::vector<size_t> &ndims)
/* returns pair for gbx index of neighbour in the
{backwards, forwards} y direction given a gridbox with
gbxidx='idx' in a cartesian domain. Treatment of neighbours
for gridboxes at the edges of the domain is either finite
(null neighbour) or periodic (cyclic neighbour) */
{
  const auto nznx = (unsigned int)ndims.at(0) * ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
  // return finitedomain_nghbrs(idx, nznx, ndims.at(2));
  return periodicdomain_nghbrs(idx, nznx, ndims.at(2));
}

KOKKOS_INLINE_FUNCTION double
coord3_beyondz(const double coord3,
               const double lim1,
               const double lim2)
/* return value is new coord for a superdroplet given that
coord3 exceedes the domain's lower or upper boundary in z direction.
(ie. coord3 is below the lower edge of the lowest gridboxes
in the z direction, or coord3 is above the upper edge of highest
gridboxes in the z direction) */
{
  return coordbeyond_finitedomain(coord3, lim1, lim2);
  // return coordbeyond_periodicdomain(coord3, lim1, lim2);
};

KOKKOS_INLINE_FUNCTION double
coord1_beyondx(const double coord1,
               const double lim1,
               const double lim2)
/* return value is new coord for a superdroplet given
that coord1 exceedes the domain's backwardsmost boundary
in x direction, or given that coord1 exceedes the
domain's forwardmost boundary in x direction */
{
  // return coordbeyond_finitedomain(coord1, lim1, lim2);
  return coordbeyond_periodicdomain(coord1, lim1, lim2);
};

KOKKOS_INLINE_FUNCTION double
coord2_beyondy(const double coord2,
               const double lim1,
               const double lim2)
/* return value is new coord for a superdroplet given
that coord2 exceedes the domain's edge/boundary in y
leftwards direction, or given that coord2 exceedes the
domain's edge/boundary in y rightwards direction */
{
  // return coordbeyond_finitedomain(coord2, lim1, lim2);
  return coordbeyond_periodicdomain(coord2, lim1, lim2);
};

#endif // CARTESIANBOUNDARYCONDS_HPP
