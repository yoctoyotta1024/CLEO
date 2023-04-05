// Author: Clara Bayley
// File: read_gbxboundaries.cpp
/* initialisatin of Maps4GridBoxes
struct from binary file */

#include "read_gbxboundaries.hpp"

GridBoxBoundaries read_gbxboundaries(std::string_view gridfile)
/* read metadata and data in binary file called 'gridfile', then
return GridBoxBoundaries instance created from that data */
{
  std::ifstream file(open_binary(gridfile));

  std::vector<VarMetadata> meta(metadata_from_binary(file));

  VarMetadata var(meta.at(0));
  file.seekg(var.b0, std::ios::beg);
  std::vector<unsigned int> gbxidxs(var.nvar, 0);
  binary_into_buffer<unsigned int>(file, gbxidxs);

  var = meta.at(1);
  file.seekg(var.b0, std::ios::beg);
  std::vector<double> gbxbounds(var.nvar, 0);
  binary_into_buffer<double>(file, gbxbounds);

  file.close();

  if (6*gbxidxs.size() != gbxbounds.size())
  {
    std::string errormsg = "sizes of gbxidxs and gbxbounds vectors"
                           " read from gridfile not consistent";
    throw std::invalid_argument(errormsg);
  }

  return GridBoxBoundaries{gbxidxs, gbxbounds};
}

double GridBoxBoundaries::gridboxarea(const unsigned int idx) const
/* calculates horizontal area of gridbox using boundaries
  corresponding to gridbox with gbxidx=idx. First finds position
  of first gbxbound (zmin) from position of idx in gbxidxs */
{
  const size_t pos = find_idx_in_gbxidxs(idx) * 6;

  const double deltax = gbxbounds[pos + 3] - gbxbounds[pos + 2]; // xmax - xmin
  const double deltay = gbxbounds[pos + 5] - gbxbounds[pos + 4]; // ymax - ymin

  return deltax * deltay;
}

double GridBoxBoundaries::gridboxvol(const unsigned int idx) const
/* calculates volume of gridbox using boundaries corresponding to
gridbox with gbxidx=idx. First finds position of first gbxbound (zmin)
for that gridbox from position of idx in gbxidxs */
{
  const size_t pos = find_idx_in_gbxidxs(idx) * 6;

  const double deltaz = gbxbounds[pos + 1] - gbxbounds[pos];     // zmax - zmin
  const double deltax = gbxbounds[pos + 3] - gbxbounds[pos + 2]; // xmax - xmin
  const double deltay = gbxbounds[pos + 5] - gbxbounds[pos + 4]; // ymax - ymin

  return deltaz * deltax * deltay;
}

size_t GridBoxBoundaries::find_idx_in_gbxidxs(const unsigned int idx) const
/* returns position in gbxidxs vector
where idx is found or raises error*/
{
  size_t pos = 0;
  for (auto gbxidx : gbxidxs)
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