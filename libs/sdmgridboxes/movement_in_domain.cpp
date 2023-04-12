
// Author: Clara Bayley
// File: "movement_in_domain.cpp"
/* Implementation of functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#include "./movement_in_domain.hpp"

/* ----- function called internally ----- */
void set_superdrop_gbxindex(const Maps4GridBoxes &mdlmaps,
                            const std::pair<double, double> zbounds,
                            const std::pair<double, double> xbounds,
                            const std::pair<double, double> ybounds,
                            const unsigned int gbxindex,
                            SuperdropWithGbxindex &SDinGBx);
/* For each direction, first check if gridbox index associated
with the superdrop in SDinGBx needs to change. If it does, implement
change by calling correct function for changing the sd_gbxindex to a
neighbouring gridbox's index in that direction */

char flag_tochange_sdgbxindex(const double coord,
                              const std::pair<double, double> bounds);
/* Determines value of the is_change flag used to signal if
the gridboxindex associated with a superdrop needs to change
and if so, in which direction the superdroplet needs to move.
returned char is: 'n' = no change, 'f' = move to forward neighbour,
'b' = move to backward neighbour */
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
    const std::pair<double, double>
        zbounds((*mdlmaps.idx2bounds_z.find(gbx.gbxindex)).second);

    const std::pair<double, double>
        xbounds((*mdlmaps.idx2bounds_x.find(gbx.gbxindex)).second);

    const std::pair<double, double>
        ybounds((*mdlmaps.idx2bounds_y.find(gbx.gbxindex)).second);

    for (auto &SDinGBx : gbx.span4SDsinGBx)
    {
      sdmmotion.move_superdroplet(w, u, v, SDinGBx.superdrop);
      //SDinGBx.superdrop.coord3 += 1.0;

      set_superdrop_gbxindex(mdlmaps, zbounds, xbounds, ybounds,
                             gbx.gbxindex, SDinGBx);
    }
  }
  
  exchange_superdroplets_between_gridboxes(SDsInGBxs, gridboxes);
}

void set_superdrop_gbxindex(const Maps4GridBoxes &mdlmaps,
                            const std::pair<double, double> zbounds,
                            const std::pair<double, double> xbounds,
                            const std::pair<double, double> ybounds,
                            const unsigned int gbxindex,
                            SuperdropWithGbxindex &SDinGBx)
/* For each direction, first check if gridbox index associated
with the superdrop in SDinGBx needs to change. If it does, implement
change by calling correct function for changing the sd_gbxindex to a
neighbouring gridbox's index in that direction */
{
  const char zflag(flag_tochange_sdgbxindex(SDinGBx.superdrop.coord3, zbounds));
  const char xflag(flag_tochange_sdgbxindex(SDinGBx.superdrop.coord1, xbounds));
  const char yflag(flag_tochange_sdgbxindex(SDinGBx.superdrop.coord2, ybounds));
 
  if (zflag == 'b')
  {
    SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_zdown(gbxindex);
  }
  else if (zflag == 'f')
  {
    SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_zup(gbxindex);
  }

  if (xflag == 'b')
  {
    SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_xbehind(gbxindex);
  }
  else if (xflag == 'f')
  {
    SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_xinfront(gbxindex);
  }

  if (yflag == 'b')
  {
    SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_yleft(gbxindex);
  }
  else if (yflag == 'f')
  {
    SDinGBx.sd_gbxindex = mdlmaps.get_neighbour_yright(gbxindex);
  } 
}

char flag_tochange_sdgbxindex(const double coord,
                              const std::pair<double, double> bounds)
/* Determines value of the is_change flag used to signal if
the gridboxindex associated with a superdrop needs to change
and if so, in which direction the superdroplet needs to move.
returned char is: 'n' = no change, 'f' = move to forward neighbour,
'b' = move to backward neighbour */
{
  const double llim = bounds.first;
  const double ulim = bounds.second;

  if (coord < llim)
  {
    return 'b'; // char to signal move SD backward a gridbox
  }
  else if (coord >= ulim)
  {
    return 'f'; // char to signal move SD forward a gridbox
  }
  else
  {
    return 'n'; // char to signal no change to SD gridbox
  }
}