// Author: Clara Bayley
// File: superdrops_with_gridboxes.hpp
/* Header file for functions involved in handling
the SuperdropWithGbxindex instances (see superdrop.hpp
for definition of this struct). Four Functions can be
called externally: 1) for creating a vector of
these objects by reading a binary file containing
intial superdrop (SD) data, 2) for sorting the
vector based on a value in each struct, 3) for changing 
an sd_gbxindex to neighbouring gridboxes' gbxindex, 
4) for printing member variables of a SDinGBx instance */

#ifndef SUPERDROPS_WITH_GRIDBOXES_HPP
#define SUPERDROPS_WITH_GRIDBOXES_HPP

#include <vector>
#include <string>
#include <string_view>
#include <memory>
#include <map>
#include <utility>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "./maps4gridboxes.hpp"
#include "claras_SDconstants.hpp"
#include "initialisation/read_initsuperdrops.hpp"
#include "superdrop_solver/superdrop.hpp"

namespace dlc = dimless_constants;

std::vector<SuperdropWithGbxindex>
superdrops_from_initSDsfile(std::string_view initSDs_filename,
                            const int nSDsvec,
                            const int SDnspace,
                            const std::shared_ptr<const SoluteProperties> solute,
                            const Maps4GridBoxes &mdlmaps);
/* reads initsuperdrop file for superdroplets' initial properties. Uses this data
to create 'nSDsvec' no. of SuperdropletWithGridbox instances in a vector
where all the superdroplets have the same solute properties, "solute".
Uses the coordinates of each superdroplet to set the value of the sd_gbxindex
associated with each superdroplet in the SuperdropletWithGridbox struct */

void sdgbxindex_to_neighbour(const Maps4GridBoxes &mdlmaps,
                                SuperdropWithGbxindex &SDinGBx);
/* first check if gridbox index associated with the superdrop
in SDinGBx needs to change. If it does, implement change by
calling correct function for changing the sd_gbxindex to a
neighbouring gridbox's index in a particular direction.
The direction is given by the value of the is_change flag */

inline void sort_superdrops_via_gridboxindex(std::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* uses the value of sd_gbxindex within each SuperdropWithGbxindex
struct to sort the vector from lowest sd_gbxindex to highest. Sorting
of objects with same value of sd_gbxindex can take any order */
{
  auto compare = [](SuperdropWithGbxindex &a, SuperdropWithGbxindex &b)
  {
    return (a.sd_gbxindex) < (b.sd_gbxindex);
  };

  std::sort(SDsInGBxs.begin(), SDsInGBxs.end(), compare);
}

inline void print_SDinGBx(const SuperdropWithGbxindex SDinGBx)
{
  std::cout << "SD " << SDinGBx.superdrop.id.value
              << ": " << SDinGBx.sd_gbxindex << ", " << SDinGBx.superdrop.eps
              << ", " << SDinGBx.superdrop.radius << ", " << SDinGBx.superdrop.m_sol
              << ", " << SDinGBx.superdrop.coord3 << ", " << SDinGBx.superdrop.coord1
              << ", " << SDinGBx.superdrop.coord2 << "\n";
}

#endif // SUPERDROPS_WITH_GRIDBOXES_HPP