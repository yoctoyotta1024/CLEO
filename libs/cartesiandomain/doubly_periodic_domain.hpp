/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: doubly_periodic_domain.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Definition of the Domain Boundary Conditions to use for Cartesian GridBox Maps, Motion of
 * Super-Droplets and MoveSupersInDomain
 */

#ifndef LIBS_CARTESIANDOMAIN_DOUBLY_PERIODIC_DOMAIN_HPP_
#define LIBS_CARTESIANDOMAIN_DOUBLY_PERIODIC_DOMAIN_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <vector>

#include "cartesiandomain/domainboundaries.hpp"

/* _Note:_ Doubly Periodic Domain Boundary Conditions are defined as:
  z: FINITE   (see cartesian_coord3nghbrs & boundarycond_coord3)
  x: PERIODIC (see cartesian_coord1nghbrs & boundarycond_coord1)
  y: PERIODIC (see cartesian_coord2nghbrs & boundarycond_coord2)
*/
namespace DoublyPeriodicDomain {

/* returns pair for gbx index of neighbour in the
{backwards, forwards} z direction given a gridbox with
gbxidx='idx' in a cartesian domain. Treatment of neighbours
for gridboxes at the edges of the domain is either finite
(null neighbour) or periodic (cyclic neighbour) */
inline Kokkos::pair<unsigned int, unsigned int> cartesian_coord3nghbrs(
    const unsigned int idx, const std::vector<size_t> &ndims) {
  return finitedomain_nghbrs(idx, 1, ndims.at(0));
  // return periodicdomain_nghbrs(idx, 1, ndims.at(0));
}

/* returns pair for gbx index of neighbour in the
{backwards, forwards} x direction given a gridbox with
gbxidx='idx' in a cartesian domain. Treatment of neighbours
for gridboxes at the edges of the domain is either finite
(null neighbour) or periodic (cyclic neighbour) */
inline Kokkos::pair<unsigned int, unsigned int> cartesian_coord1nghbrs(
    const unsigned int idx, const std::vector<size_t> &ndims) {
  const auto nz = (unsigned int)ndims.at(0);  // no. gridboxes in z direction
  // return finitedomain_nghbrs(idx, nz, ndims.at(1));
  return periodicdomain_nghbrs(idx, nz, ndims.at(1));
}

/* returns pair for gbx index of neighbour in the
{backwards, forwards} y direction given a gridbox with
gbxidx='idx' in a cartesian domain. Treatment of neighbours
for gridboxes at the edges of the domain is either finite
(null neighbour) or periodic (cyclic neighbour) */
inline Kokkos::pair<unsigned int, unsigned int> cartesian_coord2nghbrs(
    const unsigned int idx, const std::vector<size_t> &ndims) {
  const auto nznx = (unsigned int)ndims.at(0) *
                    ndims.at(1);  // no. gridboxes in z direction * no. gridboxes in x direction
  // return finitedomain_nghbrs(idx, nznx, ndims.at(2));
  return periodicdomain_nghbrs(idx, nznx, ndims.at(2));
}

/* return value is new coord for a superdroplet given that
coord3 exceedes the domain's lower or upper boundary in z direction.
(ie. coord3 is below the lower edge of the lowest gridboxes
in the z direction, or coord3 is above the upper edge of highest
gridboxes in the z direction) */
KOKKOS_INLINE_FUNCTION double boundarycond_coord3(const double coord3, const double lim1,
                                                  const double lim2) {
  return coordbeyond_finitedomain(coord3, lim1, lim2);
  // return coordbeyond_periodicdomain(coord3, lim1, lim2);
}

/* return value is new coord for a superdroplet given
that coord1 exceedes the domain's backwardsmost boundary
in x direction, or given that coord1 exceedes the
domain's forwardmost boundary in x direction */
KOKKOS_INLINE_FUNCTION double boundarycond_coord1(const double coord1, const double lim1,
                                                  const double lim2) {
  // return coordbeyond_finitedomain(coord1, lim1, lim2);
  return coordbeyond_periodicdomain(coord1, lim1, lim2);
}

/* return value is new coord for a superdroplet given
that coord2 exceedes the domain's edge/boundary in y
leftwards direction, or given that coord2 exceedes the
domain's edge/boundary in y rightwards direction */
KOKKOS_INLINE_FUNCTION double boundarycond_coord2(const double coord2, const double lim1,
                                                  const double lim2) {
  // return coordbeyond_finitedomain(coord2, lim1, lim2);
  return coordbeyond_periodicdomain(coord2, lim1, lim2);
}

}  // namespace DoublyPeriodicDomain

#endif  // LIBS_CARTESIANDOMAIN_DOUBLY_PERIODIC_DOMAIN_HPP_
