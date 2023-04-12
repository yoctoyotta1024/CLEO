
// Author: Clara Bayley
// File: "superdropsmotion.cpp"
/* Implementation of functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#include "./superdropsmotion.hpp"

void sdmmotion(const Maps4GridBoxes &mdlmaps,
               std::vector<SuperdropWithGbxindex> &SDsInGBxs,
               std::vector<GridBox> &gridboxes)
/* update superdroplet positions according to their motion and
move superdroplets between gridboxes by changing their associated
gridboxindex if necessary, then (re)sorting SDsInGBxs vector and
updating spans4SDsInGbx for each gridbox */
{
  exchange_superdroplets_between_gridboxes(mdlmaps, SDsInGBxs, gridboxes);
}

void exchange_superdroplets_between_gridboxes(const Maps4GridBoxes &mdlmaps,
                                              std::vector<SuperdropWithGbxindex> &SDsInGBxs,
                                              std::vector<GridBox> &gridboxes)
/* move superdroplets between gridboxes by changing their associated
gridboxindex if necessary, then (re)sorting SDsInGBxs vector and
updating spans4SDsInGbx for each gridbox */
{
  change_superdroplets_gridboxindex(mdlmaps, gridboxes);

  sort_superdrops_via_gridboxindex(SDsInGBxs);

  set_gridboxes_superdropletspan(gridboxes, SDsInGBxs);

  // for (auto gbx: gridboxes)
  // {
  //   gbx.iscorrect_span_for_gbxindex(mdlmaps);
  // }
}

void change_superdroplets_gridboxindex(const Maps4GridBoxes &mdlmaps,
                                       std::vector<GridBox> &gridboxes)
/* first check if superdrop's associated gridboxindex (sd_gbxindex)
needs to change. If it does, implement change by calling correct
function for changing the sd_gbxindex to a neighbouring gridbox's index
in a particular direction. The direction is given by the value of
the is_change flag */
{
  for (auto &gbx : gridboxes)
  {
    for (auto &SDinGBx : gbx.span4SDsinGBx)
    {
      sdgbxindex_to_neighbour(mdlmaps, SDinGBx); // see superdrops_with_gridboxes.cpp
    }
  }
}