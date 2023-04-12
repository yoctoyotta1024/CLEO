// Author: Clara Bayley
// File: "movement_in_domain.hpp"
/* Header file for functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#ifndef MOVEMENT_IN_DOMAIN_HPP
#define MOVEMENT_IN_DOMAIN_HPP

#include <vector>
#include <span>

#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"
#include "./superdrops_in_gridboxes.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/sdmmotion.hpp"

void move_superdrops_in_domain(const Maps4GridBoxes &mdlmaps,
                               const SdmMotion &sdmmotion,
                               std::vector<SuperdropWithGbxindex> &SDsInGBxs,
                               std::vector<GridBox> &gridboxes);
/* Move superdroplets that are in gridboxes including exchange
between gridboxes if necessary. First update superdroplet positions
according to their motion and then move superdroplets between
gridboxes by changing their associated gridboxindex as appropriate.
Final step is (re)sorting SDsInGBxs vector and updating
spans4SDsInGbx for each gridbox */

/* ----- function called internally ----- */
void exchange_superdroplets_between_gridboxes(const Maps4GridBoxes &mdlmaps,
                                              std::vector<SuperdropWithGbxindex> &SDsInGBxs,
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
/* ------------------------------------------------------ */

#endif // MOVEMENT_IN_DOMAIN_HPP