// Author: Clara Bayley
// File: "maps4gridboxes.cpp"
/* functionality for creating and using
the map between a grid box index
and its coordinate boundaries */

#include "maps4gridboxes.hpp"

std::pair<double, double> numeric_limit_bounds()
{
  return {-std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};
}

Maps4GridBoxes::Maps4GridBoxes(const unsigned int SDnspace,
                                   std::string_view gridfile)
/* initilaises idx2bounds_[i] maps (for i = x, y or z) which map
from every gridbox index to its boundaries in domain coordinates.
Also initialises idx2vol map whose values are the volume of a gridbox
given the gridbox's index as key. The keys of idx2bounds_[i] map's
are also gridbox indexes. The corresponding value is that gridbox's
{upper boundary, lower boundary}. In a non-3D case, coordinates of the
gridbox boundaries for unused dimensions are the min/max possible
doubles of computer (numerical limits), however the volume remains
finite. E.g. In the 0-D case, the idx2bounds maps have 1 {key, value}
for gridbox 0 which are the upper and lower numerical limits,
whilst the volume is determind by reading the gridfile */
{
  const GridBoxBoundaries gfb(read_gbxboundaries(gridfile, SDnspace));

  gbxidxs = gfb.gbxidxs;  
  maxidx = *std::max_element(gfb.gbxidxs.begin(), gfb.gbxidxs.end());

  if (SDnspace == 0)
  {
    const double domainvol = get_0Ddomainvol_from_gridfile(gfb);
    set_0Dmodel_maps(domainvol);
  }

  else if (SDnspace == 1)
  {
    set_1Dmodel_maps(gfb);
  }

  else if (SDnspace == 2)
  {
    set_2Dmodel_maps(gfb);
  }

  else if (SDnspace == 3)
  {
    set_3Dmodel_maps(gfb);
  }

  else
  {
    const std::string errmsg("SDnspace > 3, no method exists "
                             " for constructing Maps4GridBoxes object");
    throw std::invalid_argument(errmsg);
  }
}

void Maps4GridBoxes::set_0Dmodel_maps(const double domainvol)
/* set idx2bounds_[i] maps to numeical limits. Set volume
 map using coords read from gridfile */
{
  ndims = {1, 1, 1};

  idx2bounds_z[0] = numeric_limit_bounds();
  idx2bounds_x[0] = numeric_limit_bounds();
  idx2bounds_y[0] = numeric_limit_bounds();
  
  idx2vol[0] = domainvol; // dimensionless volume of 0D model

  idx2nghbour_z[0] = {0, 0}; // 'periodic' BCs in non-existent dimensions 
  idx2nghbour_x[0] = {0, 0};
  idx2nghbour_y[0] = {0, 0};
}

void Maps4GridBoxes::set_1Dmodel_maps(const GridBoxBoundaries &gfb)
/* Set idx2bounds x and y maps to numerical limits. Set z and volume
maps using coords from gridfile. It is assumed that for a gridbox with
it's index at position 'p' in the gfb.gbxidxs vector, the
[zmin, zmax] coords of that gridbox are at [pos, pos+1] in the
gfb.gbxidxs vector, where pos = p*6 */
{
  ndims = {gfb.ndims.at(0), 1, 1};
   
  size_t pos = 0;
  for(auto idx : gfb.gbxidxs)
  {
    idx2bounds_x[idx] = numeric_limit_bounds();
    idx2bounds_y[idx] = numeric_limit_bounds();

    const double zlow = gfb.gbxbounds[pos];
    const double zup = gfb.gbxbounds[pos+1];
    idx2bounds_z[idx] = {zlow, zup};

    const double vol = (zup - zlow) * gfb.gridboxarea(idx);
    idx2vol[idx] = vol;

    idx2nghbour_z[idx] = znghbours_cartesian(idx, gfb.gbxidxs); 
    idx2nghbour_x[idx] = {idx, idx}; // 'periodic' BCs in non-existent dimensions
    idx2nghbour_y[idx] = {idx, idx};
    
    pos += 6;
  }
}

void Maps4GridBoxes::set_2Dmodel_maps(const GridBoxBoundaries &gfb)
/* Set idx2bounds y map to numerical limits. Set z, x and volume
maps using coords from gridfile. It is assumed that for a gridbox
with it's index at position 'p' in the gfb.gbxidxs
vector, the [zmin, zmax, xmin, xmax] coords of that gridbox are
at [pos, pos+1, pos+2, pos+3] in the gfb.gbxidxs
vector, where pos = p*6 */
{
  ndims = {gfb.ndims.at(0), gfb.ndims.at(1), 1};

  size_t pos = 0;
  for(auto idx : gfb.gbxidxs)
  {
    idx2bounds_y[idx] = numeric_limit_bounds(); 

    const double zlow = gfb.gbxbounds[pos];
    const double zup = gfb.gbxbounds[pos+1];
    idx2bounds_z[idx] = {zlow, zup};

    const double xlow = gfb.gbxbounds[pos+2];
    const double xup = gfb.gbxbounds[pos+3];
    idx2bounds_x[idx] = {xlow, xup};

    const double deltay = gfb.gbxbounds[pos+5] - gfb.gbxbounds[pos+4]; 
    const double vol = (zup - zlow) * (xup - xlow) * deltay;
    idx2vol[idx] = vol;

    idx2nghbour_z[idx] = znghbours_cartesian(idx, gbxidxs); 
    idx2nghbour_x[idx] = xnghbours_cartesian(idx, gbxidxs);
    idx2nghbour_y[idx] = {idx, idx}; // 'periodic' BCs in non-existent dimensions
    
    pos += 6;
  }
}

void Maps4GridBoxes::set_3Dmodel_maps(const GridBoxBoundaries &gfb)
/*Set z, x, y and volume maps using coords from gridfile. It is assumed
that for a gridbox with it's index at position 'p' in the
gfb.gbxidxs vector, the [zmin, zmax, xmin, xmax, ymin, ymax]
coords of that gridbox are at [pos, pos+1, pos+2, pos+3, pos+4, pos+5]
in the gfb.gbxidxs vector, where pos = p*6 */
{
  ndims = {gfb.ndims.at(0), gfb.ndims.at(1), gfb.ndims.at(2)};
  
  size_t pos = 0;
  for(auto idx : gfb.gbxidxs)
  {
    const double zlow = gfb.gbxbounds[pos];
    const double zup = gfb.gbxbounds[pos+1];
    idx2bounds_z[idx] = {zlow, zup};

    const double xlow = gfb.gbxbounds[pos+2];
    const double xup = gfb.gbxbounds[pos+3];
    idx2bounds_x[idx] = {xlow, xup};

    const double ylow = gfb.gbxbounds[pos+4];
    const double yup = gfb.gbxbounds[pos+5];
    idx2bounds_y[idx] = {ylow, yup};

    const double vol = (zup - zlow) * (xup - xlow) * (yup - ylow);
    idx2vol[idx] = vol;

    idx2nghbour_z[idx] = znghbours_cartesian(idx, gbxidxs); 
    idx2nghbour_x[idx] = xnghbours_cartesian(idx, gbxidxs);
    idx2nghbour_y[idx] = ynghbours_cartesian(idx, gbxidxs);
    
    pos += 6;
  }
}

std::pair<unsigned int,
          unsigned int>
Maps4GridBoxes::znghbours_cartesian(const unsigned int idx,
                                     const std::vector<
                                         unsigned int> &gbxidxs)
/* returns gbx indexes of {upwards, downwards} neighbour
of gridbox with index idx. End points return 
max unsigned int value. */
{
  unsigned int zup = idx+1;
  if (zup > maxidx)
  {
    zup = -1; // no neighbour above gbx with largest idx
  }

  const unsigned int zdown = std::max(-1, (int)idx - 1); // no neighbour below gbx with lowest idx

  return {zup, zdown};
}

std::pair<unsigned int,
          unsigned int>
Maps4GridBoxes::xnghbours_cartesian(const unsigned int idx,
                    const std::vector<unsigned int> &gbxidxs)
{

  const unsigned int nz = ndims.at(0); // no. gridboxes in z direction

  return { idx + nz, idx - nz};
}

std::pair<unsigned int,
          unsigned int>
Maps4GridBoxes::ynghbours_cartesian(const unsigned int idx,
                    const std::vector<unsigned int> &gbxidxs)
{
  const unsigned int nznx = ndims.at(0) * ndims.at(1); // no. gridboxes in z direction * no. gridboxes in x direction
  return { idx + nznx, idx - nznx};
}