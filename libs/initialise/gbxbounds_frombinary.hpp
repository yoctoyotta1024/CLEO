/*
 * ----- CLEO -----
 * File: gbxbounds_frombinary.hpp
 * Project: initialise
 * Created Date: Wednesday 1st November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 2nd November 2023
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

public:
  std::vector<size_t> ndims;         // number of gridboxes in [coord3, coord1, coord2] dimensions
  std::vector<unsigned int> gbxidxs; // gridbox indexes
  std::vector<double> gbxbounds;     // corresponding [coord3 {l, u}, coord1 {l, u}, coord2 {l, u}] lower and upper coordinate boundaries

  GbxBoundsFromBinary(const unsigned int nspacedims,
                      std::string_view grid_filename);

  double gbxarea_fromgridfile(const unsigned int idx) const;
  /* calculates horizontal (x-y planar) area of gridbox using boundaries
   corresponding to gridbox with gbxidx=idx. First finds position
   of first gbxbound (zmin) from position of idx in gbxidxs */

  double gbxvol_fromgridfile(const unsigned int idx) const;
  /* calculates volume of gridbox using boundaries corresponding to
  gridbox with gbxidx=idx. First finds position of first gbxbound (zmin)
  for that gridbox from position of idx in gbxidxs */                  
};

#endif // GBXBOUNDS_FROMBINARY_HPP
