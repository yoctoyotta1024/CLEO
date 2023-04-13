
// Author: Clara Bayley
// File: "movement_in_domain.cpp"
/* Implementation of functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#include "./movement_in_domain.hpp"

/* ----- function called internally ----- */
unsigned int update_superdrop_gbxindex(const Maps4GridBoxes &mdlmaps,
                                       const unsigned int gbxindex,
                                       const std::pair<double, double> zbounds,
                                       const std::pair<double, double> xbounds,
                                       const std::pair<double, double> ybounds,
                                       const Superdrop &superdrop);
/* For each direction, first check if gridbox index associated
with the superdrop in SDinGBx needs to change. If it does, implement
change by calling correct function for changing the sd_gbxindex to a
neighbouring gridbox's index in that direction */
/* -------------------------------------- */

void move_superdrops_in_domain(const Maps4GridBoxes &mdlmaps,
                               const SdmMotion &sdmmotion,
                               std::vector<SuperdropWithGbxindex> &SDsInGBxs,
                               std::vector<GridBox> &gridboxes)
/* Move superdroplets that are in gridboxes including exchange
between gridboxes if necessary. First update superdroplet positions
according to their motion and then move superdroplets between
gridboxes by changing their associated gridboxindex as appropriate.
Final step is (re)sorting SDsInGBxs vector and updating
spans4SDsInGbx for each gridbox */
{
  const double w = 0.0; // TODO: get winds from gridbox state
  const double u = 0.0;
  const double v = 0.0;

  for (auto &gbx : gridboxes)
  {
    const auto zbounds(mdlmaps.get_bounds_z(gbx.gbxindex));
    const auto xbounds(mdlmaps.get_bounds_x(gbx.gbxindex));
    const auto ybounds(mdlmaps.get_bounds_y(gbx.gbxindex));

    for (auto &SDinGBx : gbx.span4SDsinGBx)
    {
      sdmmotion.move_superdroplet(w, u, v, SDinGBx.superdrop);

      SDinGBx.sd_gbxindex = update_superdrop_gbxindex(mdlmaps, gbx.gbxindex,
                                                      zbounds, xbounds, ybounds,
                                                      SDinGBx.superdrop);
    }
  }

  exchange_superdroplets_between_gridboxes(SDsInGBxs, gridboxes);
}

unsigned int update_superdrop_gbxindex(const Maps4GridBoxes &mdlmaps,
                                       const unsigned int gbxindex,
                                       const std::pair<double, double> zbounds,
                                       const std::pair<double, double> xbounds,
                                       const std::pair<double, double> ybounds,
                                       const Superdrop &superdrop)
/* For each direction (z, then x, then y), mdlmaps's forward and backward 
get_neighbour functions are passed into changeindex_ifcoord_outofbounds
along with superdroplet's coord and the gridbox bounds for that direction.
If coord not within bounds, changeindex_ifcoord_outofbounds is used to
return a new value of sd_gbxindex via calling the appropriate get_neighbour
function. After algorithm for z, then x, then y directions are complete,
resultant sd_gbxindex is returned. */
{
  auto bindfunc = [](const Maps4GridBoxes mdlmaps, const auto funcref)
  {
    return std::bind(funcref, mdlmaps, std::placeholders::_1); // cast sd_gbxindex to *signed* int
  };

  unsigned int sd_gbxindex(gbxindex);

  const auto zdown = bindfunc(mdlmaps, &Maps4GridBoxes::get_neighbour_zdown);
  const auto zup = bindfunc(mdlmaps, &Maps4GridBoxes::get_neighbour_zup);
  sd_gbxindex = changeindex_ifcoord_outofbounds(sd_gbxindex, zbounds,
                                                   superdrop.coord3,
                                                   zdown, zup);

  const auto xbehind = bindfunc(mdlmaps, &Maps4GridBoxes::get_neighbour_xbehind);
  const auto xinfront = bindfunc(mdlmaps, &Maps4GridBoxes::get_neighbour_xinfront);
  sd_gbxindex = changeindex_ifcoord_outofbounds(sd_gbxindex, xbounds,
                                                   superdrop.coord1,
                                                   xbehind, xinfront);

  const auto yleft = bindfunc(mdlmaps, &Maps4GridBoxes::get_neighbour_yleft);
  const auto yright = bindfunc(mdlmaps, &Maps4GridBoxes::get_neighbour_yright);
  sd_gbxindex = changeindex_ifcoord_outofbounds(sd_gbxindex, ybounds,
                                                   superdrop.coord2,
                                                   yleft, yright);

  return sd_gbxindex;
}