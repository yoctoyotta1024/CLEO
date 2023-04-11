// Author: Clara Bayley
// File: "maps4gridboxes.hpp"
/* Header file for functions
related to creating and using
maps between a gridbox indicies
and domain coordinates */

/* their coordinate boundaries
or between a girdbox and
the index of its neighbour in a
given direction */

#ifndef MAPS4GRIDBOXES_HPP
#define MAPS4GRIDBOXES_HPP

#include <map>
#include <utility>
#include <string_view>
#include <string>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <array>

#include "initialisation/readbinary.hpp"
#include "initialisation/read_gbxboundaries.hpp"

struct Maps4GridBoxes
{
private:
  void set_0Dmodel_maps(const double domainvol);
  /* set vol map using coords read from gridfile */

  void set_1Dmodel_maps(const GridBoxBoundaries &gfb);
  /* set z and vol maps using coords from gridfile */

  void set_2Dmodel_maps(const GridBoxBoundaries &gfb);
  /* Set z, x and volume maps using coords from gridfile */

  void set_3Dmodel_maps(const GridBoxBoundaries &gfb);
  /* Set z, x, y and volume maps using coords from gridfile. */

public:
  std::vector<unsigned int> gbxidxs;                              // vector of all gridbox indexes in domain
  std::map<unsigned int, std::pair<double, double>> idx2bounds_z; // coord limits to each gridbox given its index
  std::map<unsigned int, std::pair<double, double>> idx2bounds_x;
  std::map<unsigned int, std::pair<double, double>> idx2bounds_y;
  std::map<unsigned int, double> idx2vol; // volume of gridbox given its index

  std::map<unsigned int, std::pair<unsigned int, unsigned int>> idx2nghbour_z; // neigbouring gbxindex to each gridbox given its gbxindex
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> idx2nghbour_x;
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> idx2nghbour_y;

  Maps4GridBoxes(const unsigned int SDnspace, std::string_view gridfile);
  /* initilaises idx2bounds_[i] maps (for i = x, y or z) which map
  from every gridbox index to its boundaries in domain coordinates.
  Also initialises idx2vol map whose values are the volume of a gridbox
  given the gridbox's index as key. The keys of idx2bounds_[i] map's
  are also gridbox indexes. The corresponding value is that gridbox's
  {upper boundary, lower boundary}. In a non-3D case, coordinates of the
  gridbox boundaries for unused dimensions are the min/max possible
  doubles of computer (numerical limits), however the volume remains
  finite. E.g. In the 0-D case, the idx2bounds maps have 1 {key, value}
  for gridbox 0 which are the upper and lower numerical limits,
  whilst the volume is determind by reading the gridfile */

  inline unsigned int get_neighbour_zup(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards z, ie. upwards direction */
  {
    return (*idx2nghbour_z.find(gbxindex)).second.first;
  }

  inline unsigned int get_neighbour_zdown(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards z, ie. downwards direction */
  {
    return (*idx2nghbour_z.find(gbxindex)).second.second;
  }

  inline unsigned int get_neighbour_yright(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards y direction, ie. right */
  {
    return (*idx2nghbour_y.find(gbxindex)).second.first;
  }

  inline unsigned int get_neighbour_yleft(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the backwards y direction, ie. left */
  {
    return (*idx2nghbour_y.find(gbxindex)).second.second;
  }

  inline unsigned int get_neighbour_xforward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards x direction, ie. out of page */
  {
    return (*idx2nghbour_x.find(gbxindex)).second.first;
  }

  inline unsigned int get_neighbour_xbackward(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the backwards x direction, ie. into page */
  {
    return (*idx2nghbour_x.find(gbxindex)).second.second;
  }
};

struct CartesianNeighbourIndexes
{
private:
  unsigned int maxidx;                     // largest value gridbox index
  std::array<size_t, 3> ndims = {0, 0, 0}; // number of gridboxes in [z,x,y] directions

  std::pair<unsigned int, unsigned int>
  handle_finitedomain_nghbours(const unsigned int forward,
                               const unsigned int backward) const;
  /* retunrs {forward, backward} gridbox neighbours with
  treatment of neighbours as if bounds of domain are finite.
  Means that no neighbour exists above/below highest/lowest gbxindex.
  For non-existent neighbours, max unsigned int value is returned,
  ie. neighbour backwards for gridboxes with backward<0 is maximum unsigned int,
  while neighbour forwards for gridboxes with forward>maxidx
  is maximm unsigned int */

  std::pair<unsigned int, unsigned int>
  handle_periodicdomain_nghbours(const unsigned int forward,
                                 const unsigned int backward) const;                             

public:
  CartesianNeighbourIndexes(const unsigned int maxidx,
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
    return handle_finitedomain_nghbours(idx + 1, idx - 1);
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
    return handle_finitedomain_nghbours(idx + nz, idx - nz);
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
    return handle_finitedomain_nghbours(idx + nznx, idx - nznx);
  }
};

#endif // MAPS4GRIDBOXES_HPP