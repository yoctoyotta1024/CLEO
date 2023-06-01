// Author: Clara Bayley
// File: "gridbox.cpp"
/* Functionality relevant to
a gridbox */

#include "gridbox.hpp"

KOKKOS_FUNCTION
GridBox::GridBox(const unsigned int ii,
                 const Maps4GridBoxes &gbxmaps,
                 const InstallDetectors &dtrs,
                 Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs)
    : gbxindex(ii),
      state(gbxmaps.get_volume(gbxindex)),
      detectors(dtrs(gbxindex, gbxmaps))
/* Volume in Thermostate set using Map4GridBoxes
idx2vol map (via get_volume function). Other ThermoState variables
are default behaviour initialised. */
{
  // print_statevolume();
  
  set_span(SDsInGBxs);
  iscorrect_span_for_gbxindex(gbxmaps);
}

KOKKOS_FUNCTION void GridBox::print_statevolume()
/* print's dimensionless value for gridbox state's 
volume. Also prints true volume = volume * COORD0^3 [m^3] */
{
  const double vol = state.get_volume();
  std::cout << "dimensionless volume = " << vol
            << "\nie. VOLUME = "
            << vol * pow(dlc::COORD0, 3.0) << "m^3\n";
}

KOKKOS_FUNCTION
void GridBox::set_span(Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs)
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

KOKKOS_FUNCTION
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

KOKKOS_FUNCTION
void GridBox::iscoord_within_bounds(const std::pair<double, double> bounds,
                                    const double coord)
{
  const double llim = bounds.first;
  const double ulim = bounds.second;

  if (coord < llim || coord >= ulim)
  {
    const std::string err("superdrop coord: "+std::to_string(coord)+
                              " lies outside its gridbox's bounds ["+
                              std::to_string(llim)+", "+
                              std::to_string(ulim)+"]");
    throw std::invalid_argument(err);
  }
}

KOKKOS_FUNCTION Kokkos::vector<GridBox>
create_gridboxes(const Maps4GridBoxes &gbxmaps,
                 Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs)
/* create domain as a vector of grid boxes such that each grid box
is initialised with a labels from gbxmaps.gbxidxs, and a span of the
superdroplet 'SDsInGbxs', and an (uninitialised) thermodynamic state. */
{ 
  sort_superdrops_via_gridboxindex(SDsInGBxs);
  
  Kokkos::vector<GridBox> gridboxes;
  for (auto ii : gbxmaps.gbxidxs)
  {
    gridboxes.push_back(GridBox(ii, gbxmaps, SDsInGBxs));
  }

  return gridboxes;
}

KOKKOS_FUNCTION
void set_superdroplets_to_wetradius(Kokkos::vector<GridBox> &gridboxes)
/* for each gridbox, set the radius of each superdroplet (SD) to
whichever is larger out of their dry radius or equlibrium wet radius
(given the relative humidity (s_ratio) and temperature of the gridbox).
If relh > maxrelh = 0.95, set each SD's radius to their
equilibrium radius at relh = maxrelh = 0.95 */
{

  const double maxrelh = 0.95;

  for (auto &gbx : gridboxes)
  {
    const double temp = gbx.state.temp;
    const double psat = saturation_pressure(temp);
    const double supersat = supersaturation_ratio(gbx.state.press,
                                                  gbx.state.qvap, psat);
    const double s_ratio = std::min(maxrelh, supersat);
    
    for (auto &SDinGBx : gbx.span4SDsinGBx)
    {
      const double wetr(SDinGBx.superdrop.equilibrium_wetradius(s_ratio,
                                                                temp));
      const double dryr(SDinGBx.superdrop.get_dry_radius());
      SDinGBx.superdrop.radius = std::max(dryr, wetr);
    }
  }
}