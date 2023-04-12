// Author: Clara Bayley
// File: superdrops_in_gridboxes.hpp
/* Header file for functions involved in handling
vector of SuperdropWithGbxindex instances (see superdrop.hpp
for definition of this struct) associated with Gridboxes
defined by relations of gbxindex in Maps4GridBoxes.
Note: some hidden functions called internally are
defined in .cpp implementaiton file */

#ifndef SUPERDROPS_IN_GRIDBOXES_HPP
#define SUPERDROPS_IN_GRIDBOXES_HPP

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

inline void print_SDinGBx(const SuperdropWithGbxindex SDinGBx)
{
  std::cout << "SD " << SDinGBx.superdrop.id.value
              << ": " << SDinGBx.sd_gbxindex << ", " << SDinGBx.superdrop.eps
              << ", " << SDinGBx.superdrop.radius << ", " << SDinGBx.superdrop.m_sol
              << ", " << SDinGBx.superdrop.coord3 << ", " << SDinGBx.superdrop.coord1
              << ", " << SDinGBx.superdrop.coord2 << "\n";
}



#endif // SUPERDROPS_IN_GRIDBOXES_HPP