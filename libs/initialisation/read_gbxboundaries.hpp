// Author: Clara Bayley
// File: read_gbxboundaries.hpp
/* Header for initialisation of Maps4GridBoxes
struct from binary file */

#ifndef READ_GBXBOUNDARIES_HPP
#define READ_GBXBOUNDARIES_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string_view>

#include "./readbinary.hpp"

struct GridBoxBoundaries
/* holds vectors containing gridbox indicies and the corresponding
coords of the [zmin, zmax, zmin, xmax, ymin, ymax] boundaries of 
that gridbox which are read from gridfile and used in
construction of Maps4GridBoxes */
{
  std::vector<unsigned int> gbxidxs; // gridbox indicies
  std::vector<double> gbxbounds; // corresponding (z,x,y) coords of max and min boundaries
  
  // double domainarea() const
  // /* returns horizontal area of entire domain */
  // {
  //   const double xdelta = std::abs(xhalf.front() - xhalf.back());
  //   const double ydelta = std::abs(yhalf.front() - yhalf.back());

  //   return xdelta * ydelta;
  // }

  double gridboxvol(const unsigned int gbxidx) const
  /* calculates volume of gridbox using boundaries corresponding to
  gridbox with gbxidx=idx */
  {
    const double zdelta = std::abs(zhalf.front() - zhalf.back());

    return zdelta * domainarea();
  }
};

GridBoxBoundaries read_gbxboundaries(std::string_view gridfile);
/* read metadata and data in binary file called 'gridfile', then
return GridBoxBoundaries instance created from that data */

inline double get_0Ddomainvol_from_gridfile(std::string_view gridfile)
/* return the volume of the 0th gridbox by reading the 'gridfile'
binary. This is the domian volume in the 0D (1 gridbox) model */
{
  const GridBoxBoundaries gbxbounds(read_gbxboundaries(gridfile));

  return gbxbounds.gridboxvol(0);
}

#endif // READ_GBXBOUNDARIES_HPP 