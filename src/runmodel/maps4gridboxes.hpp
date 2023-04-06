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
  std::map<unsigned int, std::pair<double, double>> idx2bounds_z; // coord limits to each gridbox given its index
  std::map<unsigned int, std::pair<double, double>> idx2bounds_x;
  std::map<unsigned int, std::pair<double, double>> idx2bounds_y;
  std::map<unsigned int, double> idx2vol; // volume of gridbox given its index

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

  inline unsigned int get_gridboxneighbour_up(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in upwards direction */
  {
    return gbxindex + 1;
  }

  inline unsigned int get_gridboxneighbour_down(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in downwards direction */
  {
    return gbxindex - 1;
  }
};

#endif // MAPS4GRIDBOXES_HPP