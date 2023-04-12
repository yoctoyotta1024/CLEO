
// Author: Clara Bayley
// File: "movement_in_domain.cpp"
/* Implementation of functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#include "./movement_in_domain.hpp"

/* ----- function called internally ----- */
void sdgbxindex_to_neighbour(const Maps4GridBoxes &mdlmaps,
                                SuperdropWithGbxindex &SDinGBx);
/* first check if gridbox index associated with the superdrop
in SDinGBx needs to change. If it does, implement change by
calling correct function for changing the sd_gbxindex to a
neighbouring gridbox's index in a particular direction.
The direction is given by the value of the is_change flag */

int flag_tochange_sdgbxindex(const SuperdropWithGbxindex &SDinGBx,
                             const std::map<unsigned int,
                                            std::pair<double, double>> &idx2bounds_z);
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
  const double w = 0.0; // TODO: get winds from gridbox
  const double u = 0.0;
  const double v = 0.0;
  
  for (auto &gbx : gridboxes)
  {
    for (auto &SDinGBx : gbx.span4SDsinGBx)
    {
      sdmmotion.move_superdroplet(w, u, v, SDinGBx.superdrop);

      sdgbxindex_to_neighbour(mdlmaps, SDinGBx);
    }
  }
  
  exchange_superdroplets_between_gridboxes(SDsInGBxs, gridboxes);
}

void sdgbxindex_to_neighbour(const Maps4GridBoxes &mdlmaps,
                             SuperdropWithGbxindex &SDinGBx)
/* first check if gridbox index associated with the superdrop
in SDinGBx needs to change. If it does, implement change by
calling correct function for changing the sd_gbxindex to a
neighbouring gridbox's index in a particular direction.
The direction is given by the value of the is_change flag */
{
  enum MoveSuperdrop
  {
    No,
    Down,
    Up,
    Left,
    Right,
    Forwards,
    Backwards
  };

  int is_change = flag_tochange_sdgbxindex(SDinGBx, mdlmaps.idx2bounds_z);
  
  if (is_change != No)
  {
    if (is_change == Down)
    {
      SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_zdown(SDinGBx.sd_gbxindex);
    }
    else if (is_change == Up)
    {
      SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_zup(SDinGBx.sd_gbxindex);
    }
    else
    {
      throw std::invalid_argument("method to change SD index for"
                                  "given is_change flag is not defined");
    }
  }
}

int flag_tochange_sdgbxindex(const SuperdropWithGbxindex &SDinGBx,
                             const std::map<unsigned int, std::pair<double, double>> &idx2bounds_z)
/* Determines value of the is_change flag used to signal if
the gridboxindex associated with a superdrop needs to change
and if so, in which direction the superdroplet needs to move */
{
  const double coord = SDinGBx.superdrop.coord3;
  const double llim = (*idx2bounds_z.find(SDinGBx.sd_gbxindex)).second.first;
  const double ulim = (*idx2bounds_z.find(SDinGBx.sd_gbxindex)).second.second;

  if (coord < llim)
  {
    return 1; // enum for moving SD index down a gridbox
  }
  else if (coord >= ulim)
  {
    return 2; // enum for moving SD index up a gridbox
  }
  else
  {
    return 0; // enum for not changing SD index
  }
}