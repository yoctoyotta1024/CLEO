// Author: Clara Bayley
// File: superdropwithgbxindex.hpp
/* Header file for functions involved in handling
vector of SuperdropWithGbxindex instances (see superdrop.hpp
for definition of this struct).
Note: some hidden functions called internally are
defined in .cpp implementaiton file */

#ifndef SUPERDROPWITHGBXINDEX_HPP
#define SUPERDROPWITHGBXINDEX_HPP

#include <vector>
#include <string>
#include <string_view>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <span>

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>

#include "../claras_SDconstants.hpp"
#include "initialisation/read_initsuperdrops.hpp"
#include "superdrop_solver/superdrop.hpp"

namespace dlc = dimless_constants;

Kokkos::vector<SuperdropWithGbxindex>
create_superdrops_from_initSDsfile(std::string_view initSDs_filename,
                            const int nSDsvec,
                            const int SDnspace,
                            const std::shared_ptr<const SoluteProperties> solute);
/* reads initsuperdrop file for superdroplets' initial properties. Uses this data
to create 'nSDsvec' no. of SuperdropletWithGridbox instances in a vector
where all the superdroplets have the same solute properties, "solute".
Uses the coordinates of each superdroplet to set the value of the sd_gbxindex
associated with each superdroplet in the SuperdropletWithGridbox struct */

inline void print_SDinGBx(const SuperdropWithGbxindex &SDinGBx)
{
  std::cout << "SD " << SDinGBx.superdrop.id.value
              << ": " << SDinGBx.sd_gbxindex << ", " << SDinGBx.superdrop.eps
              << ", " << SDinGBx.superdrop.radius << ", " << SDinGBx.superdrop.m_sol
              << ", " << SDinGBx.superdrop.coord3 << ", " << SDinGBx.superdrop.coord1
              << ", " << SDinGBx.superdrop.coord2 << "\n";
}

inline void sort_superdrops_via_gridboxindex(Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs)
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

inline std::span<SuperdropWithGbxindex>
remove_outofdomain_superdrops(std::span<SuperdropWithGbxindex> &span4SDsinGBx)
/* sorts span based on sd_gbxindex and then obtains a
subspan which excludes superdroplets with OUTOFDOMAIN
value as sd_gbxindex (ie. max possible value) */
{
  /* 1. Sort span based on sd_gbxindexes */
  auto compare = [](SuperdropWithGbxindex &a, SuperdropWithGbxindex &b)
  {
    return (a.sd_gbxindex) < (b.sd_gbxindex);
  };
  std::sort(span4SDsinGBx.begin(), span4SDsinGBx.end(), compare);  

  /* 2. Find first instance where sd_gbxindex < OUTOFDOMAIN  */
  auto lowcompare = [](const SuperdropWithGbxindex &a, const unsigned int val)
  {
    return a.sd_gbxindex < val;
  };
  auto up = std::lower_bound(span4SDsinGBx.begin(), span4SDsinGBx.end(),
                             dlc::OUTOFDOMAIN, lowcompare);

  /* 3. Return subspan that is upto and excluding OUTOFDOMAIN sd_gbxindexes */  
  return {span4SDsinGBx.begin(), up};
}

#endif // SUPERDROPWITHGBXINDEX_HPP