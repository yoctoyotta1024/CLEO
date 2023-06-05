// Author: Clara Bayley
// File: "maps4gridboxes.hpp"
/* Header file for functions
related to creating and using
maps between a gridbox indicies
and domain coordinates */

#ifndef MAPS4GRIDBOXES_HPP
#define MAPS4GRIDBOXES_HPP

#include <map>
#include <utility>
#include <string_view>
#include <string>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <vector>

#include "initialisation/readbinary.hpp"
#include "initialisation/read_gbxboundaries.hpp"
#include "./cartesianneighbours.hpp"

struct Maps4GridBoxes
{
private:
  std::map<unsigned int, std::pair<double, double>> idx2bounds_z; // coord limits to each gridbox given its index
  std::map<unsigned int, std::pair<double, double>> idx2bounds_x; // value pair is {lower bound, upper bounds} for gbxindex key
  std::map<unsigned int, std::pair<double, double>> idx2bounds_y;
  std::map<unsigned int, double> idx2area; // x-y planar area of gridbox given its index
  std::map<unsigned int, double> idx2vol; // volume of gridbox given its index

  std::map<unsigned int, std::pair<unsigned int, unsigned int>> idx2nghbour_z; // neigbouring gbxindex to each gridbox given its gbxindex
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> idx2nghbour_x;
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> idx2nghbour_y;
  
  void set_0Dmodel_maps(const double domainarea, const double domainvol);
  /* set vol map using coords read from gridfile */

  void set_1Dmodel_maps(const GridBoxBoundaries &gfb);
  /* set z and vol maps using coords from gridfile */

  void set_2Dmodel_maps(const GridBoxBoundaries &gfb);
  /* Set z, x and volume maps using coords from gridfile */

  void set_3Dmodel_maps(const GridBoxBoundaries &gfb);
  /* Set z, x, y and volume maps using coords from gridfile. */

  void check_ngridboxes() const;

public:
  std::vector<unsigned int> gbxidxs;                              // vector of all gridbox indexes in domain
  std::array<size_t, 3> ndims;                                    // number of gridboxes in [z,x,y] directions
  size_t ngridboxes;

  Maps4GridBoxes(const unsigned int SDnspace, std::string_view gridfile);
  /* initilaises idx2bounds_[i] maps (for i = x, y or z) which map
  from every gridbox index to its boundaries in domain coordinates.
  Also initialises idx2area and idx2vol maps whose values are the
  area and volume of a gridbox given the gridbox's index as key.
  The keys of idx2bounds_[i] map's are also gridbox indexes. The
  corresponding value is that gridbox's {upper boundary, lower boundary}.
  In a non-3D case, coordinates of the gridbox boundaries for unused
  dimensions are the min/max possible doubles of computer (numerical
  limits), however the area and volume remain finite. E.g. In the 0-D
  case, the idx2bounds maps have 1 {key, value} for gridbox 0 which
  are the upper and lower numerical limits, whilst the volume is 
  determined by reading the gridfile */

  std::pair<double, double> get_bounds_z(const unsigned int gbxidx) const
  {
    return (*idx2bounds_z.find(gbxidx)).second;
  }

  std::pair<double, double> get_bounds_x(const unsigned int gbxidx) const
  {
    return (*idx2bounds_x.find(gbxidx)).second;
  }

  std::pair<double, double> get_bounds_y(const unsigned int gbxidx) const
  {
    return (*idx2bounds_y.find(gbxidx)).second;
  }

  double get_area(const unsigned int gbxidx) const
  {
    return (*idx2area.find(gbxidx)).second;
  }

  double get_volume(const unsigned int gbxidx) const
  {
    return (*idx2vol.find(gbxidx)).second;
  }

  unsigned int get_neighbour_zdown(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards z, ie. downwards direction */
  {
    return (*idx2nghbour_z.find(gbxindex)).second.second;
  }

  unsigned int get_neighbour_zup(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards z, ie. upwards direction */
  {
    return (*idx2nghbour_z.find(gbxindex)).second.first;
  }
  
  unsigned int get_neighbour_xbehind(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the backwards x direction, ie. into page */
  {
    return (*idx2nghbour_x.find(gbxindex)).second.second;
  }

  unsigned int get_neighbour_xinfront(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards x direction, ie. out of page */
  {
    return (*idx2nghbour_x.find(gbxindex)).second.first;
  }

  unsigned int get_neighbour_yleft(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the backwards y direction, ie. left */
  {
    return (*idx2nghbour_y.find(gbxindex)).second.second;
  }

  unsigned int get_neighbour_yright(unsigned int gbxindex) const
  /* given gridbox index, return index of neighbouring
  gridbox in the forwards y direction, ie. right */
  {
    return (*idx2nghbour_y.find(gbxindex)).second.first;
  }
};


#endif // MAPS4GRIDBOXES_HPP