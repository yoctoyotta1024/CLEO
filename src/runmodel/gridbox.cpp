// Author: Clara Bayley
// File: "gridbox.cpp"
/* Functionality relevant to
a gridbox */

#include "gridbox.hpp"

GridBox::GridBox(const unsigned int ii,
                 const Maps4GridBoxes &mdlmaps,
                 std::vector<SuperdropWithGridbox> &SDsInGBxs)
    : gbxindex(ii), state()
{
  set_statevolume(mdlmaps.idx2vol);
  
  set_span(SDsInGBxs);
  iscorrect_span_for_gbxindex(mdlmaps);
}

void GridBox::set_statevolume(const std::map<unsigned int, double> &idx2vol)
/* set dimensionless value for gridbox state's 
volume using Map4GridBoxes idx2vol map.
True volume = state.volume * COORD0^3 [m^3] */
{
  state.volume = (*idx2vol.find(gbxindex)).second;
  std::cout << "dimensionless volume = " << state.volume
            << "\nie. VOLUME = "
            << state.volume * pow(dlc::COORD0, 3.0) << "m^3\n";
}

void GridBox::set_span(std::vector<SuperdropWithGridbox> &SDsInGBxs)
/* assumes SDsInGBxs is ordered based on sd_gbxindex
from lowest to highest. Finds first and last SDWithGBx that has 
sd_gbxindex matching gbxindex in order to set span4SDsinGBx. */
{
  auto upcompare = [](const int val, const SuperdropWithGridbox &a)
  {
    return val < (int)a.sd_gbxindex; // cast sd_gbxindex to *signed* int
  };

  auto low = std::upper_bound(SDsInGBxs.begin(), SDsInGBxs.end(),
                              gbxindex - 1, upcompare);
  auto up = std::upper_bound(SDsInGBxs.begin(), SDsInGBxs.end(),
                             gbxindex, upcompare);

  span4SDsinGBx = {low, up};  
}

void GridBox::iscorrect_span_for_gbxindex(const Maps4GridBoxes &mdlmaps)
{
  for (auto &SDinGBx : span4SDsinGBx)
  {
    if(SDinGBx.sd_gbxindex != gbxindex)
    {
      const std::string err = "span4SDsinGBx incorrectly set."
                              " At least one sd_gbxindex"
                              " does not match this gridbox's index (ie. " +
                              std::to_string(SDinGBx.sd_gbxindex) +
                              " != "+std::to_string(gbxindex)+")";
      throw std::invalid_argument(err);
    }
    iscoord_within_bounds(mdlmaps.idx2bounds_z, SDinGBx.superdrop.coord3);
    iscoord_within_bounds(mdlmaps.idx2bounds_x, SDinGBx.superdrop.coord1);
    iscoord_within_bounds(mdlmaps.idx2bounds_y, SDinGBx.superdrop.coord2);
  }
}

void GridBox::iscoord_within_bounds(const std::map<unsigned int,
                                                   std::pair<double,
                                                             double>> &idx2bounds,
                                    const double coord)
{
  const double llim = (*idx2bounds.find(gbxindex)).second.first;
  const double ulim = (*idx2bounds.find(gbxindex)).second.second;

  if (coord < llim || coord >= ulim)
  {
    const std::string err = "superdrop coord: "+std::to_string(coord)+
                              " lies outside its gridbox's bounds ["+
                              std::to_string(llim)+", "+
                              std::to_string(ulim)+"]";
    throw std::invalid_argument(err);
  }
}

std::vector<GridBox> create_gridboxes(const Maps4GridBoxes &mdlmaps,
                                      std::vector<SuperdropWithGridbox> &SDsInGBxs)
/* create domain as a vector of grid boxes such that each grid box
is initialised with a labels from mdlmaps.gbxidxs, and a span of the
superdroplet 'SDsInGbxs', and an (uninitialised) thermodynamic state. */
{ 
  sort_superdrops_via_gridboxindex(SDsInGBxs);
  
  std::vector<GridBox> gridboxes;
  for (auto ii : mdlmaps.gbxidxs)
  {
    gridboxes.push_back(GridBox(ii, mdlmaps, SDsInGBxs));
  }

  return gridboxes;
}