/*
 * ----- CLEO -----
 * File: gbxbounds_frombinary.hpp
 * Project: initialise
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 3rd November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functions for reading gridbox boundaries from
 * a binary file (used to then create a
 * map from gbxindexes to gridbox boundaries
 * for CLEO SDM e.g. a CartesianMaps)
 */

#ifndef GBXBOUNDS_FROMBINARY_HPP
#define GBXBOUNDS_FROMBINARY_HPP

#include <string_view>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <array>
#include <algorithm>
#include <iterator>

#include "Kokkos_Core.hpp"
#include "Kokkos_Pair.hpp"

#include "./readbinary.hpp"

struct GbxBoundsFromBinary
/* holds vectors containing gridbox indexes and their
corresponding [zmin, zmax, zmin, xmax, ymin, ymax]
coordinate boundaries which are read from gridfile
and used in construction of GridboxMaps */
{
private:
  void is_nspacedims_compatible(const unsigned int nspacedims) const;
  
  bool check_0Dmodel_gbxbounds() const;
  bool check_1Dmodel_gbxbounds() const;
  bool check_2Dmodel_gbxbounds() const;
  bool check_3Dmodel_gbxbounds() const;

  size_t find_idx_in_gbxidxs(const unsigned int idx) const;

public:
  std::vector<size_t> ndims;         // number of gridboxes in [coord3, coord1, coord2] dimensions
  std::vector<unsigned int> gbxidxs; // gridbox indexes
  std::vector<double> gbxbounds;     // corresponding [coord3 {l, u}, coord1 {l, u}, coord2 {l, u}] lower and upper coordinate boundaries

  GbxBoundsFromBinary(const unsigned int nspacedims,
                      std::string_view grid_filename);
  
  Kokkos::pair<double, double>
  get_coord3gbxbounds(const unsigned int idx) const;
  /* returns coord3 {lower, upper} gridbox bounds 
  from position in gbxbounds vector which corresponds
  to position in gbxidxs where gbxidx = idx */

  Kokkos::pair<double, double>
  get_coord1gbxbounds(const unsigned int idx) const;
  /* returns coord1 {lower, upper} gridbox bounds 
  from position in gbxbounds vector which corresponds
  to position in gbxidxs where gbxidx = idx */

  Kokkos::pair<double, double>
  get_coord2gbxbounds(const unsigned int idx) const;
  /* returns coord2 {lower, upper} gridbox bounds 
  from position in gbxbounds vector which corresponds
  to position in gbxidxs where gbxidx = idx */

  double gbxarea_fromgridfile(const unsigned int idx) const;
  /* calculates horizontal (x-y planar) area of gridbox
  using boundaries corresponding to gridbox with gbxidx=idx. */

  double gbxvol_fromgridfile(const unsigned int idx) const;
  /* calculates volume of gridbox using boundaries
  corresponding to gridbox with gbxidx=idx. */

  size_t get_ngbxs() const
  /* returns total number of gridboxes = product of dimentions */
  {
    return ndims.at(0) * ndims.at(1) * ndims.at(2);
  }
};

#endif // GBXBOUNDS_FROMBINARY_HPP
