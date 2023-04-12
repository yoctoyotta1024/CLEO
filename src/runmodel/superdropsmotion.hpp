// Author: Clara Bayley
// File: "superdropsmotion.hpp"
/* Header file for functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#ifndef SUPERDROPSMOTION_HPP
#define SUPERDROPSMOTION_HPP

#include <vector>

#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"
#include "./superdrops_with_gridboxes.hpp"
#include "superdrop_solver/superdrop.hpp"

void sdmmotion(const Maps4GridBoxes &mdlmaps,
               std::vector<SuperdropWithGridbox> &SDsInGBxs,
               std::vector<GridBox> &gridboxes);
/* update superdroplet positions according to their motion and
move superdroplets between gridboxes by changing their associated
gridboxindex if necessary, then (re)sorting SDsInGBxs vector and
updating spans4SDsInGbx for each gridbox */

void exchange_superdroplets_between_gridboxes(const Maps4GridBoxes &mdlmaps,
                                              std::vector<SuperdropWithGridbox> &SDsInGBxs,
                                              std::vector<GridBox> &gridboxes);
/* move superdroplets between gridboxes by changing their associated
gridboxindex if necessary, then (re)sorting SDsInGBxs vector and
updating spans4SDsInGbx for each gridbox */

void change_superdroplets_gridboxindex(const Maps4GridBoxes &mdlmaps,
                                       std::vector<GridBox> &gridboxes);
/* first check if superdrop's associated gridboxindex (sd_gbxindex)
needs to change. If it does, implement change by calling correct
function for changing the sd_gbxindex to a neighbouring gridbox's index
in a particular direction. The direction is given by the value of
the is_change flag */

#endif // SUPERDROPSMOTION_HPP