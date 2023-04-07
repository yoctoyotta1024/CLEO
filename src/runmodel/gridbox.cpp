// Author: Clara Bayley
// File: "gridbox.cpp"
/* Functionality relevant to
a gridbox */

#include "gridbox.hpp"

GridBox::GridBox(const unsigned int ii,
                 const std::map<unsigned int, double> &idx2vol,
                 std::vector<SuperdropWithGridbox> &SDsInGBxs)
    : gbxindex(ii), state()
{
  set_span(SDsInGBxs);
  set_statevolume(idx2vol);
}

void GridBox::set_span(std::vector<SuperdropWithGridbox> &SDsInGBxs)
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

void GridBox::set_statevolume(const std::map<unsigned int, double> &idx2vol)
/* set dimensionless value for gridbox state's 
volume using Map4GridBoxes idx2vol map.
True volume = state.volume * COORD0^3 [m^3] */
{
  state.volume = (*idx2vol.find(gbxindex)).second;
  std::cout << "dimensionless volume = " << state.volume << "\n";
  std::cout << "ie. VOLUME = " << state.volume * pow(dlc::COORD0, 3.0) << "m^3\n";
}

std::vector<GridBox> create_gridboxes(const size_t num_gridboxes,
                                      const std::map<unsigned int, double> &idx2vol,
                                      std::vector<SuperdropWithGridbox> &SDsInGBxs)
/* create domain as a vector of grid boxes such that each grid box 
is initialised with a label (ii), a superdroplet vector with
superdroplets created from the SDinitialisation csv file, 
and an (uninitialised) thermodynamic state. */                                   
{ 
  sort_superdrops_via_gridboxindex(SDsInGBxs);
  
  std::vector<GridBox> gridboxes;
  for (unsigned int ii = 0; ii < num_gridboxes; ++ii)
  {
    gridboxes.push_back(GridBox(ii, idx2vol, SDsInGBxs));
  }

  return gridboxes;
}