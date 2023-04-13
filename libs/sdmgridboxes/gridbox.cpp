// Author: Clara Bayley
// File: "gridbox.cpp"
/* Functionality relevant to
a gridbox */

#include "gridbox.hpp"

GridBox::GridBox(const unsigned int ii,
                 const Maps4GridBoxes &gbxmaps,
                 std::vector<SuperdropWithGbxindex> &SDsInGBxs)
    : gbxindex(ii), state(gbxmaps.get_volume(gbxindex))
/* Volume in Thermostate set using Map4GridBoxes
idx2vol map (via get_volume function). Other ThermoState variables
are default behaviour initialised. */
{
  print_statevolume();
  
  set_span(SDsInGBxs);
  iscorrect_span_for_gbxindex(gbxmaps);
}

void GridBox::print_statevolume()
/* print's dimensionless value for gridbox state's 
volume. Also prints true volume = state.volume * COORD0^3 [m^3] */
{
  std::cout << "dimensionless volume = " << state.volume
            << "\nie. VOLUME = "
            << state.volume * pow(dlc::COORD0, 3.0) << "m^3\n";
}

void GridBox::set_span(std::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* assumes SDsInGBxs is ordered based on sd_gbxindex
from lowest to highest. Finds first and last SDWithGBx that has 
sd_gbxindex matching gbxindex in order to set span4SDsinGBx. */
{
  auto lowcompare = [](const SuperdropWithGbxindex &a, const unsigned int val)
  {
    return a.sd_gbxindex < val; // cast sd_gbxindex to *signed* int
  };

  auto upcompare = [](const unsigned int val, const SuperdropWithGbxindex &a)
  {
    return val < a.sd_gbxindex; // cast sd_gbxindex to *signed* int
  };

  auto low = std::lower_bound(SDsInGBxs.begin(), SDsInGBxs.end(),
                              gbxindex, lowcompare);
  auto up = std::upper_bound(SDsInGBxs.begin(), SDsInGBxs.end(),
                             gbxindex, upcompare);

  span4SDsinGBx = {low, up};  
}

void GridBox::iscorrect_span_for_gbxindex(const Maps4GridBoxes &gbxmaps)
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
    iscoord_within_bounds(gbxmaps.get_bounds_z(gbxindex), SDinGBx.superdrop.coord3);
    iscoord_within_bounds(gbxmaps.get_bounds_x(gbxindex), SDinGBx.superdrop.coord1);
    iscoord_within_bounds(gbxmaps.get_bounds_y(gbxindex), SDinGBx.superdrop.coord2);
  }
}

void GridBox::iscoord_within_bounds(const std::pair<double, double> bounds,
                                    const double coord)
{
  const double llim = bounds.first;
  const double ulim = bounds.second;

  if (coord < llim || coord >= ulim)
  {
    const std::string err = "superdrop coord: "+std::to_string(coord)+
                              " lies outside its gridbox's bounds ["+
                              std::to_string(llim)+", "+
                              std::to_string(ulim)+"]";
    throw std::invalid_argument(err);
  }
}

std::vector<GridBox> create_gridboxes(const Maps4GridBoxes &gbxmaps,
                                      std::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* create domain as a vector of grid boxes such that each grid box
is initialised with a labels from gbxmaps.gbxidxs, and a span of the
superdroplet 'SDsInGbxs', and an (uninitialised) thermodynamic state. */
{ 
  sort_superdrops_via_gridboxindex(SDsInGBxs);
  
  std::vector<GridBox> gridboxes;
  for (auto ii : gbxmaps.gbxidxs)
  {
    gridboxes.push_back(GridBox(ii, gbxmaps, SDsInGBxs));
  }

  return gridboxes;
}