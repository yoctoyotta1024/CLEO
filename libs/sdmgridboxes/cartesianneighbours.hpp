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

#include <utility>
#include <array>
#include <vector>
#include <limits>

struct CartesianNeighbourGBxIndexes
{
private:
  unsigned int maxidx;                     // largest value gridbox index
  std::array<size_t, 3> ndims = {0, 0, 0}; // number of gridboxes in [z,x,y] directions

  std::pair<unsigned int, unsigned int>
  handle_finitedomain_nghbours(const unsigned int idx,
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
  handle_periodicdomain_nghbours(const unsigned int idx,
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
  'handle_XXX_nghbours' function */
  {
    return handle_finitedomain_nghbours(idx, 1, ndims.at(0));
  }

  std::pair<unsigned int, unsigned int>
  xnghbours_cartesian(const unsigned int idx,
                      const std::vector<unsigned int> &gbxidxs) const
  /* returns pair of gbx indexes for {infront, behind} neighbour
  of a gridbox with index 'idx'. Treatment of neighbours for
  gridboxes at edges of domain is determined by the
  'handle_XXX_nghbours' function */
  {
    const unsigned int nz = ndims.at(0); // no. gridboxes in z direction
    return handle_finitedomain_nghbours(idx, nz, ndims.at(1));
  }

  std::pair<unsigned int, unsigned int>
  ynghbours_cartesian(const unsigned int idx,
                      const std::vector<unsigned int> &gbxidxs) const
  /* returns pair of gbx indexes for {right, left} neighbour
  of a gridbox with index 'idx'. Treatment of neighbours for
  gridboxes at edges of domain is determined by the
  'handle_XXX_nghbours' function */
  {
    const unsigned int nznx = ndims.at(0) * ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
    return handle_finitedomain_nghbours(idx, nznx, ndims.at(2));
  }
};

#endif // CARTESIANNEIGHBOURS_HPP