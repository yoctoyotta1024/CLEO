/*
 * ----- CLEO -----
 * File: gbxbounds_frombinary.cpp
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


#include "./gbxbounds_frombinary.hpp"

GbxBoundsFromBinary::
    GbxBoundsFromBinary(const unsigned int nspacedims,
                            std::string_view grid_filename)
/* read metadata and data in binary file called 'gridfile', then
return GbxBoundsFromBinary instance created from that data */
{
  /* open file and read in the metatdata
  for all the variables in gridfile */
  std::ifstream file(open_binary(grid_filename));
  std::vector<VarMetadata> meta(metadata_from_binary(file));

  ndims = vector_from_binary<size_t>(file, meta.at(0));
  gbxidxs = vector_from_binary<unsigned int>(file, meta.at(1)); 
  gbxbounds = vector_from_binary<double>(file, meta.at(2)); 

  file.close();

  if (6 * gbxidxs.size() != gbxbounds.size() && gbxbounds.size() >=6)
  {
    std::string errormsg = "sizes of gbxidxs and gbxbounds vectors"
                           " read from gridfile not consistent";
    throw std::invalid_argument(errormsg);
  }

  is_nspacedims_compatible(nspacedims);
}

void GbxBoundsFromBinary::
    is_nspacedims_compatible(const unsigned int nspacedims) const
/* check that nspacedims is consistent with ndims and then
calls appropropriate function to check if gbxbounds is also.
finally throws error if either proves inconsistent */                                       
{
  bool isgood = false;

  if (nspacedims == 0)
  {
    isgood = check_0Dmodel_gbxbounds();
  }

  else if (nspacedims == 1 &&
           ndims.at(1) == 1 && ndims.at(2) == 1)
  {
    isgood = check_1Dmodel_gbxbounds();
  }

  else if (nspacedims == 2 && ndims.at(2) == 1 )
  {
    // 2D model should have constant y coords
    isgood = check_2Dmodel_gbxbounds();
  }

  else if (nspacedims == 3)
  {
    isgood = check_3Dmodel_gbxbounds(); 
  }

  else
  {
    std::string err("ndims from gridfile and/or SDnspace not valid");
    throw std::invalid_argument(err);
  }

  if (isgood == false)
  {
    std::string err = "gridbox bounds read from gridfile "
                      "not compatible with nspacedims = " +
                      std::to_string(nspacedims);
    throw std::invalid_argument(err);
  }
}

bool GbxBoundsFromBinary::check_0Dmodel_gbxbounds() const
/* returns true if data for gridbox boundaries, gbxbounds,
is compatible with 0-D model. Criteria is that 0-D model
has 1 gridbox and hence 6 values in gbxbounds */
{
  if (gbxbounds.size() == 6 &&
      ndims.at(0) == 1 && ndims.at(1) == 1 &&
      ndims.at(2) == 1)
  {
    return true;
  }

  return false;
}

bool GbxBoundsFromBinary::check_1Dmodel_gbxbounds() const
/* returns true if data for gridbox boundaries, gbxbounds,
is compatible with 1D model. Criteria is that x and y
coords of all gridbox boundaries are the same. */  
{
  const std::array<double, 4> bounds0{gbxbounds[2],
                                      gbxbounds[3],
                                      gbxbounds[4],
                                      gbxbounds[5]};

  size_t ii(2);                      // start at 0th gridbox's xlow (i.e. 3rd value in gbxbounds)
  while (ii + 4 <= gbxbounds.size()) // loop over each gridbox's bounds
  {
    for (int j = 0; j < 4; ++j)
    {
      // check x and y bounds are same as given by bounds0
      if (bounds0[j] != gbxbounds[ii + j])
      {
        return false;
      }
    }

    ii += 6;
  }

  return true;
}

bool GbxBoundsFromBinary::check_2Dmodel_gbxbounds() const
/* returns true if data for gridbox boundaries,
gbxbounds, is compatible with 2D model. Criteria is
that y coords of all gridbox boundaries are the same. */  
{
  const std::array<double, 2> bounds0{gbxbounds[4], gbxbounds[5]};

  size_t ii(4);                      // start at 0th gridbox's ylow (i.e. 5th value in gbxbounds)
  while (ii + 2 <= gbxbounds.size()) // loop over each gridbox's bounds
  {
    for (int j = 0; j < 2; ++j)
    {
      // check y bounds are same as given by bounds0
      if (bounds0[j] != gbxbounds[ii + j])
      {
        return false;
      }
    }

    ii += 6;
  }

  return true;
}

bool GbxBoundsFromBinary::check_3Dmodel_gbxbounds() const
/* returns true if data for gridbox boundaries,
gbxbounds, is compatible with 0-D model. Criteria
is that 3-D model should have at least 1 gridbox */
{

  if (gbxbounds.size() >= 6)
  {
    return true;
  }
  return false;
}

size_t GbxBoundsFromBinary::
    find_idx_in_gbxidxs(const unsigned int idx) const
/* returns distance (number of hops) from start of
gbxidxs vector to position where gbxidx matches idx
(ie. *it = idx) */
{
  auto it(std::find(gbxidxs.begin(), gbxidxs.end(), idx)); //iterator to idx
  size_t pos(std::distance(gbxidxs.begin(), it)); // distance from start of gbxidxs to idx

  if (pos > (gbxidxs.size()-1))
  {
    /* if pos is larger than the largest valid position
    in gbxidxs, idx has not been found so throw an error*/
    throw std::invalid_argument("idx not found in gbxidxs");
  }

  return pos;
}

Kokkos::pair<double, double>
GbxBoundsFromBinary::get_coord3gbxbounds(const unsigned int idx) const
/* returns coord3 {lower, upper} gridbox bounds
from position in gbxbounds vector which corresponds
to position in gbxidxs where gbxidx = idx */
{
  const unsigned int pos = find_idx_in_gbxidxs(idx) * 6; // position of zmin for gbxidx = idx

  return {gbxbounds[pos], gbxbounds[pos + 1]};
}

Kokkos::pair<double, double>
GbxBoundsFromBinary::get_coord1gbxbounds(const unsigned int idx) const
/* returns coord1 {lower, upper} gridbox bounds
from position in gbxbounds vector which corresponds
to position in gbxidxs where gbxidx = idx
'pos' is position of first bound (ie. zmin) for gridbox
assuming order is [zmin, zmax, xmin, xmax, ymin, ymax] */
{
  const unsigned int pos = find_idx_in_gbxidxs(idx) * 6; // position of zmin for gbxidx = idx

  return {gbxbounds[pos], gbxbounds[pos + 1]};
}

Kokkos::pair<double, double>
GbxBoundsFromBinary::get_coord2gbxbounds(const unsigned int idx) const
/* returns coord2 {lower, upper} gridbox bounds
from position in gbxbounds vector which corresponds
to position in gbxidxs where gbxidx = idx
'pos' is position of first bound (ie. zmin) for gridbox
assuming order is [zmin, zmax, xmin, xmax, ymin, ymax] */
{
  const unsigned int pos = find_idx_in_gbxidxs(idx) * 6; // position of zmin for gbxidx = idx

  return {gbxbounds[pos + 2], gbxbounds[pos + 3]};
}

Kokkos::pair<double, double>
GbxBoundsFromBinary::get_coord2gbxbounds(const unsigned int idx) const
/* returns coord2 {lower, upper} gridbox bounds
from position in gbxbounds vector which corresponds
to position in gbxidxs where gbxidx = idx. 
'pos' is position of first bound (ie. zmin) for gridbox
assuming order is [zmin, zmax, xmin, xmax, ymin, ymax] */
{
  const unsigned int pos = find_idx_in_gbxidxs(idx) * 6; // position of zmin for gbxidx = idx

  return {gbxbounds[pos + 4], gbxbounds[pos + 5]};
}

double GbxBoundsFromBinary::
    gbxarea_fromgridfile(const unsigned int idx) const
/* calculates horizontal (x-y planar) area of gridbox
using boundaries corresponding to gridbox with gbxidx=idx. */
{
  const auto xbounds(get_coord1gbxbounds(idx));
  const auto ybounds(get_coord2gbxbounds(idx));

  const double deltax(xbounds.second - xbounds.first);     // xmax - xmin
  const double deltay(ybounds.second - ybounds.first);     // ymax - ymin

  return deltax * deltay;
}

double GbxBoundsFromBinary::
    gbxvol_fromgridfile(const unsigned int idx) const
/* calculates volume of gridbox using boundaries
corresponding to gridbox with gbxidx=idx. */
{
  const auto zbounds(get_coord3gbxbounds(idx));
  const auto xbounds(get_coord1gbxbounds(idx));
  const auto ybounds(get_coord2gbxbounds(idx));

  const double deltaz(zbounds.second - zbounds.first);     // zmax - zmin
  const double deltax(xbounds.second - xbounds.first);     // xmax - xmin
  const double deltay(ybounds.second - ybounds.first);     // ymax - ymin

  return deltaz * deltax * deltay;
}

