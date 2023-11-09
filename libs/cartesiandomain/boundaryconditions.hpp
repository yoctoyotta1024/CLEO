/*
 * ----- CLEO -----
 * File: boundaryconditions.hpp
 * Project: cartesiandomain
 * Created Date: Thursday 9th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 9th November 2023
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

#ifndef BOUNDARYCONDITIONS_HPP
#define BOUNDARYCONDITIONS_HPP

/* NOTE: boundary conditions of domain are defined as:
  z: FINITE   (see cartesian_znghbrs & coord3_beyondz)
  x: PERIODIC (see cartesian_xnghbrs & coord1_beyondx)
  y: PERIODIC (see cartesian_ynghbrs & coord2_beyondy)
*/

#include <Kokkos_Core.hpp>

KOKKOS_INLNE_FUNCTION bool
at_cartesiandomainboundary(const unsigned int gbxindex,
                           const unsigned int increment,
                           const size_t ndim)
/* returns true if gbxindex for gridbox is at a boundary
of cartesian domain, given neighbouring indexes are
+- increment from gbxindex and the number of gridboxes
making up the domain in that direction (ndim) */
{
  return (gbxindex / increment) % ndim == 0;
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

#endif // BOUNDARYCONDITIONS_HPP 
