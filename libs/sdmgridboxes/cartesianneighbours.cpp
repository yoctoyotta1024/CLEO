// Author: Clara Bayley
// File: "cartesianneighbours.hpp"
/* Some functions involved in
 calculating neighbouring
gridbox indexes when initialising 
Maps4Gridboxes, and functions for
updating superdroplet coords given
peiodic boundary conditions */

#include "cartesianneighbours.hpp"

std::pair<unsigned int, unsigned int>
CartesianNeighbourGBxIndexes::finitedomain_nghbours(const unsigned int idx,
                                                        const unsigned int increment,
                                                        const unsigned int ndim) const
/* returns {forward, backward} gridbox neighbours with
treatment of neighbours as if bounds of domain are finite. This
means that no neighbour exists above/below highest/lowest gridboxes 
in a given direction. For non-existent neighbours, max unsigned int
value is returned, ie. in a given direction, index of neighbour
backwards and/or forwards of gridboxes at edge of domain is
maximum unsigned int */
{
  unsigned int forward = idx + increment;
  unsigned int backward = idx - increment;

  if ((idx/increment) % ndim == 0) // at lower edge of domain
  {
    backward = std::numeric_limits<unsigned int>::max();
  }

  if ((forward/increment) % ndim == 0) // at upper edge of domain
  {
    forward = std::numeric_limits<unsigned int>::max();
  }

  return {forward, backward};
}

std::pair<unsigned int, unsigned int>
CartesianNeighbourGBxIndexes::periodicdomain_nghbours(const unsigned int idx,
                              const unsigned int increment,
                              const unsigned int ndim) const
/* returns {forward, backward} gridbox neighbours with
treatment of neighbours as if bounds of domain are periodic. This
means that highest/lowest gridboxes in a given direction are 
neighbours. I.e.  index of neighbour forwards of gridboxes at
the uppermost edge of domain in a given direction, are the
lowermost gridboxes in that direction (and vice versa). */
{
  unsigned int forward = idx + increment;
  unsigned int backward = idx - increment;

  if ((idx/increment) % ndim == 0) // at lower edge of domain
  {
    backward = idx + (ndim-1) * increment;
  }

  if ((forward/increment) % ndim == 0) // at upper edge of domain
  {
    forward = idx - (ndim-1) * increment;
  }

  return {forward, backward};
}