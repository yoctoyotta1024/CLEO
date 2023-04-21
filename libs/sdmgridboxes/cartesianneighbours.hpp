// Author: Clara Bayley
// File: "cartesianneighbours.hpp"
/* Header file for functions
related to calculating neighbouring
gridbox indexes when initialising 
Maps4Gridboxes. Also functions updating
superdroplet coords given peiodic 
boundary conditions */

#ifndef CARTESIANNEIGHBOURS_HPP
#define CARTESIANNEIGHBOURS_HPP

#include <cstddef>
#include <vector>
#include <utility>
#include <limits>
#include <array>

/* boundary conditions of domain are defined as:
  z: FINITE   (see znghbours_cartesian, coord3_beyondzdown and coord3_beyondzup)
  x: PERIODIC (see xnghbours_cartesian, coord1_beyondxbehind and coord1_beyondxinfront)
  y: PERIODIC (see ynghbours_cartesian, coord2_beyondyleft and coord2_beyondyright)
*/

inline bool at_domainboundary(const unsigned int idx,
                              const unsigned int increment,
                               const unsigned int ndim)
/* returns true if idx for gridbox is at a domain boundary, given
neighbouring indexes are +- increment from idx and the number of
gridboxes making up the domain in that direction (ndim) */
{
  return (idx/increment) % ndim == 0;
}

struct CartesianNeighbourGBxIndexes
{
private:
  unsigned int maxidx;                     // largest value gridbox index
  std::array<size_t, 3> ndims = {0, 0, 0}; // number of gridboxes in [z,x,y] directions

  std::pair<unsigned int, unsigned int>
  finitedomain_nghbours(const unsigned int idx,
                              const unsigned int increment,
                              const unsigned int ndim) const;
  /* returns {forward, backward} gridbox neighbours with
  treatment of neighbours as if bounds of domain are finite. This
  means that no neighbour exists above/below highest/lowest gridboxes 
  in a given direction. For non-existent neighbours, max unsigned int
  value is returned, ie. in a given direction, index of neighbour
  backwards and/or forwards of gridboxes at edge of domain is
  maximum unsigned int */

  std::pair<unsigned int, unsigned int>
  periodicdomain_nghbours(const unsigned int idx,
                              const unsigned int increment,
                              const unsigned int ndim) const;                       
  /* returns {forward, backward} gridbox neighbours with
  treatment of neighbours as if bounds of domain are periodic. This
  means that highest/lowest gridboxes in a given direction are 
  neighbours. I.e.  index of neighbour forwards of gridboxes at
  the uppermost edge of domain in a given direction, are the
  lowermost gridboxes in that direction (and vice versa). */

public:
  CartesianNeighbourGBxIndexes(const unsigned int maxidx,
                            const std::array<size_t, 3> ndims)
      : maxidx(maxidx), ndims(ndims) {}

  std::pair<unsigned int, unsigned int>
  znghbours_cartesian(const unsigned int idx,
                      const std::vector<
                          unsigned int> &gbxidxs) const
  /* returns pair of gbx indexes for {upwards, downwards} neighbour
  of a gridbox with index 'idx'. Treatment of neighbours for
  gridboxes at edges of domain is determined by the
  'XXXdomain_nghbours' function */
  {
    return finitedomain_nghbours(idx, 1, ndims.at(0));
    // return periodicdomain_nghbours(idx, 1, ndims.at(0));
  }

  std::pair<unsigned int, unsigned int>
  xnghbours_cartesian(const unsigned int idx,
                      const std::vector<unsigned int> &gbxidxs) const
  /* returns pair of gbx indexes for {infront, behind} neighbour
  of a gridbox with index 'idx'. Treatment of neighbours for
  gridboxes at edges of domain is determined by the
  'XXXdomain_nghbours' function */
  {
    const unsigned int nz = ndims.at(0); // no. gridboxes in z direction
    // return finitedomain_nghbours(idx, nz, ndims.at(1));
    return periodicdomain_nghbours(idx, nz, ndims.at(1));
  }

  std::pair<unsigned int, unsigned int>
  ynghbours_cartesian(const unsigned int idx,
                      const std::vector<unsigned int> &gbxidxs) const
  /* returns pair of gbx indexes for {right, left} neighbour
  of a gridbox with index 'idx'. Treatment of neighbours for
  gridboxes at edges of domain is determined by the
  'XXXdomain_nghbours' function */
  {
    const unsigned int nznx = ndims.at(0) * ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
    // return finitedomain_nghbours(idx, nznx, ndims.at(2));
    return periodicdomain_nghbours(idx, nznx, ndims.at(2));
  }
};

inline double coordbeyond_finitedomain(const double coord)
/* Finite domain boundaries don't change superdroplet coord */
{
  return coord; // finite domain therefore don't change coord
}

inline double coord3_beyondzdown(const double coord3)
/* return value is new coord for a superdroplet given that
coord3 exceedes the domain's lower boundary in z direction
(ie. coord3 is below the minimum gridbox bound in z direction */
{
  return coordbeyond_finitedomain(coord3);
};

inline double coord3_beyondzup(const double coord3)
/* return value is new coord for a superdroplet given that
coord3 exceedes the domain's upper boundary in z direction
(ie. coord3 is above the maximum gridbox bound in z direction */
{
  return coordbeyond_finitedomain(coord3);
};


#endif // CARTESIANNEIGHBOURS_HPP