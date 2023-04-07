// Author: Clara Bayley
// File: "maps4gridboxes.cpp"
/* functionality for creating and using
the map between a grid box index
and its coordinate boundaries */

#include "maps4gridboxes.hpp"

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
  idx2bounds_z[0] = {-std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};

  idx2bounds_x[0] = {-std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};

  idx2bounds_y[0] = {-std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};

  idx2vol[0] = domainvol; // dimensionless volume of 0D model
}

void Maps4GridBoxes::set_1Dmodel_maps(const GridBoxBoundaries &gfb)
/* Set idx2bounds x and y maps to numerical limits. Set z and volume
maps using coords from gridfile. It is assumed that for a gridbox with
it's index at position 'p' in the gfb.gbxidxs vector, the
[zmin, zmax] coords of that gridbox are at [pos, pos+1] in the
gfb.gbxidxs vector, where pos = p*6 */
{
  idx2bounds_x[0] = {-std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};

  idx2bounds_y[0] = {-std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};

  size_t pos = 0;
  for(auto idx : gfb.gbxidxs)
  {
    const double zlow = gfb.gbxbounds[pos];
    const double zup = gfb.gbxbounds[pos+1];
    idx2bounds_z[idx] = {zlow, zup};

    const double vol = (zup - zlow) * gfb.gridboxarea(idx);
    idx2vol[idx] = vol;

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
  idx2bounds_y[0] = {-std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};

  size_t pos = 0;
  for(auto idx : gfb.gbxidxs)
  {
    const double zlow = gfb.gbxbounds[pos];
    const double zup = gfb.gbxbounds[pos+1];
    idx2bounds_z[idx] = {zlow, zup};

    const double xlow = gfb.gbxbounds[pos+2];
    const double xup = gfb.gbxbounds[pos+3];
    idx2bounds_x[idx] = {xlow, xup};

    const double deltay = gfb.gbxbounds[pos+5] - gfb.gbxbounds[pos+4]; 
    const double vol = (zup - zlow) * (xup - xlow) * deltay;
    idx2vol[idx] = vol;

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

    pos += 6;
  }
}