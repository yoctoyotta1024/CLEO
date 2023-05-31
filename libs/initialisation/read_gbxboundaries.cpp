// Author: Clara Bayley
// File: read_gbxboundaries.cpp
/* initialisatin of Maps4GridBoxes
struct from binary file */

#include "read_gbxboundaries.hpp"

void is_gridbounds_SDnspace_compatible(const unsigned int SDnspace,
                                       const std::vector<double> &gbxbounds,
                                       std::vector<size_t> &ndims);
/* check that data for gridbox boundaries read from gridfile is
compatible with SDnspace from config file. Throw error if not. */                                       

bool check_1Dmodel_gridbounds(const std::vector<double> &gbxbounds);
/* returns true if data for gridbox boundaries, gbxbounds,
is compatible with 1D model. Criteria is that x and y coords of
all gridbox boundaries are the same. */    

bool check_2Dmodel_gridbounds(const std::vector<double> &gbxbounds);
/* returns true if data for gridbox boundaries, gbxbounds,
is compatible with 2D model. Criteria is that y coords of
all gridbox boundaries are the same. */  

GridBoxBoundaries read_gbxboundaries(std::string_view gridfile,
                                     const unsigned int SDnspace)
/* read metadata and data in binary file called 'gridfile', then
return GridBoxBoundaries instance created from that data */
{
  /* open file and read in the metatdata
  for all the variables in gridfile */
  std::ifstream file(open_binary(gridfile));
  std::vector<VarMetadata> meta(metadata_from_binary(file));

  std::vector<size_t>
      ndims(vector_from_binary<size_t>(file, meta.at(0)));

  std::vector<unsigned int>
      gbxidxs(vector_from_binary<unsigned int>(file, meta.at(1))); 

  std::vector<double>
      gbxbounds(vector_from_binary<double>(file, meta.at(2))); 

  file.close();

  if (6 * gbxidxs.size() != gbxbounds.size() && gbxbounds.size() >=6)
  {
    std::string errormsg = "sizes of gbxidxs and gbxbounds vectors"
                           " read from gridfile not consistent";
    throw std::invalid_argument(errormsg);
  }

  is_gridbounds_SDnspace_compatible(SDnspace, gbxbounds, ndims);

  return GridBoxBoundaries{ndims, gbxidxs, gbxbounds};
}

double GridBoxBoundaries::gridboxarea(const unsigned int idx) const
/* calculates horizontal area of gridbox using boundaries
  corresponding to gridbox with gbxidx=idx. First finds position
  of first gbxbound (zmin) from position of idx in gbxidxs */
{
  const unsigned int pos = find_idx_in_gbxidxs(idx) * 6;

  const double deltax = gbxbounds[pos + 3] - gbxbounds[pos + 2]; // xmax - xmin
  const double deltay = gbxbounds[pos + 5] - gbxbounds[pos + 4]; // ymax - ymin

  return deltax * deltay;
}

double GridBoxBoundaries::gridboxvol(const unsigned int idx) const
/* calculates volume of gridbox using boundaries corresponding to
gridbox with gbxidx=idx. First finds position of first gbxbound (zmin)
for that gridbox from position of idx in gbxidxs */
{
  const unsigned int pos = find_idx_in_gbxidxs(idx) * 6;

  const double deltaz = gbxbounds[pos + 1] - gbxbounds[pos];     // zmax - zmin
  const double deltax = gbxbounds[pos + 3] - gbxbounds[pos + 2]; // xmax - xmin
  const double deltay = gbxbounds[pos + 5] - gbxbounds[pos + 4]; // ymax - ymin

  return deltaz * deltax * deltay;
}

unsigned int GridBoxBoundaries::find_idx_in_gbxidxs(const unsigned int idx) const
/* returns position in gbxidxs vector
where idx is found or raises error */
{
  unsigned int pos = 0;
  for (const auto gbxidx : gbxidxs)
  {
    if (gbxidx == idx)
    {
      break;
    }
    ++pos;
  }

  if (pos >= gbxidxs.size())
  {
    std::string errormsg = "index of gridbox, " +
                           std::to_string(idx) +
                           ", not found in gbxidxs vector";
    throw std::invalid_argument(errormsg);
  }

  return pos;
}

void is_gridbounds_SDnspace_compatible(const unsigned int SDnspace,
                                       const std::vector<double> &gbxbounds,
                                       std::vector<size_t> &ndims)
{
  bool isgood = false;

  if (SDnspace == 0)
  {
    // 0D model should have 1 gridbox, hence 6 values in gbxbounds
    if (gbxbounds.size() == 6 &&
        ndims.at(0) == 1 && ndims.at(1) == 1 &&
        ndims.at(2) == 1)
    {
      isgood = true;
    }
  }

  else if (SDnspace == 1 &&
           ndims.at(1) == 1 && ndims.at(2) == 1)
  {
    isgood = check_1Dmodel_gridbounds(gbxbounds);
  }

  else if (SDnspace == 2 && ndims.at(2) == 1 )
  {
    // 2D model should have constant y coords
    isgood = check_2Dmodel_gridbounds(gbxbounds);
  }

  else if (SDnspace == 3)
  {
    // 3D model should have at least 1 gridbox
    if (gbxbounds.size() >= 6)
    {
      isgood = true;
    }
  }

  else
  {
    throw std::invalid_argument("ndims from gridfile and/or SDnspace not valid");
  }

  if (isgood != true)
  {
    std::string errormsg = "gridbounds read from gridfile not compatible"
                           " with SDnspace = " +
                           std::to_string(SDnspace);
    throw std::invalid_argument(errormsg);
  }
}

bool check_1Dmodel_gridbounds(const std::vector<double> &gbxbounds)
/* returns true if data for gridbox boundaries, gbxbounds,
is compatible with 1D model. Criteria is that x and y coords of
all gridbox boundaries are the same. */  
{
  const std::array<double, 4> bounds0{gbxbounds[2], gbxbounds[3],
                                gbxbounds[4], gbxbounds[5]};

  size_t ii(2);                      // start at 0th gridbox's xlow (i.e. 3rd value in gbxbounds)
  while (ii + 4 <= gbxbounds.size()) // loop over each gridbox's bounds
  {
    for (int j = 0; j < 4; ++j)
    {
      // check x and y bounds are same as bounds0
      if (bounds0[j] != gbxbounds[ii + j])
      {
        return false;
      }
    }

    ii += 6;
  }

  return true;
}

bool check_2Dmodel_gridbounds(const std::vector<double> &gbxbounds)
/* returns true if data for gridbox boundaries, gbxbounds,
is compatible with 2D model. Criteria is that y coords of
all gridbox boundaries are the same. */  
{
  const std::array<double, 2> bounds0{gbxbounds[4], gbxbounds[5]};

  size_t ii(4);                      // start at 0th gridbox's ylow (i.e. 5th value in gbxbounds)
  while (ii + 2 <= gbxbounds.size()) // loop over each gridbox's bounds
  {
    for (int j = 0; j < 2; ++j)
    {
      // check y bounds are same as bounds0
      if (bounds0[j] != gbxbounds[ii + j])
      {
        return false;
      }
    }

    ii += 6;
  }

  return true;
}