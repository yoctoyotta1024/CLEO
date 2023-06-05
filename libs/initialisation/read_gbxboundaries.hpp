// Author: Clara Bayley
// File: read_gbxboundaries.hpp
/* Header for initialisation of Maps4GridBoxes
struct from binary file */

#ifndef READ_GBXBOUNDARIES_HPP
#define READ_GBXBOUNDARIES_HPP

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <string_view>
#include <stdexcept>

#include "./readbinary.hpp"

struct GridBoxBoundaries
/* holds vectors containing gridbox indicies and the corresponding
coords of the [zmin, zmax, zmin, xmax, ymin, ymax] boundaries of
that gridbox which are read from gridfile and used in
construction of Maps4GridBoxes */
{
private:
  unsigned int find_idx_in_gbxidxs(const unsigned int idx) const;

public:
  std::vector<size_t> ndims; // number of gridboxes in (z,x,y) dimensions
  std::vector<unsigned int> gbxidxs; // gridbox indicies
  std::vector<double> gbxbounds;     // corresponding (z,x,y) coords of max and min boundaries

  double gridboxarea(const unsigned int idx) const;
  /* calculates horizontal area of gridbox using boundaries
   corresponding to gridbox with gbxidx=idx. First finds position
   of first gbxbound (zmin) from position of idx in gbxidxs */

  double gridboxvol(const unsigned int idx) const;
  /* calculates volume of gridbox using boundaries corresponding to
  gridbox with gbxidx=idx. First finds position of first gbxbound (zmin)
  for that gridbox from position of idx in gbxidxs */
};

GridBoxBoundaries read_gbxboundaries(std::string_view gridfile,
                                     const unsigned int SDnspace);
/* read metadata and data in binary file called 'gridfile', then
return GridBoxBoundaries instance created from that data */

inline double get_0Ddomainarea_from_gridfile(const GridBoxBoundaries &gfb)
/* return the volume of the 0th gridbox by reading the 'gridfile'
binary. This is the domian volume in the 0D (1 gridbox) model */
{
  return gfb.gridboxarea(0);
}

inline double get_0Ddomainvol_from_gridfile(const GridBoxBoundaries &gfb)
/* return the volume of the 0th gridbox by reading the 'gridfile'
binary. This is the domian volume in the 0D (1 gridbox) model */
{
  return gfb.gridboxvol(0);
}

#endif // READ_GBXBOUNDARIES_HPP