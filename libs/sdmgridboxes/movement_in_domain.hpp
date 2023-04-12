// Author: Clara Bayley
// File: "movement_in_domain.hpp"
/* Header file for functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#ifndef MOVEMENT_IN_DOMAIN_HPP
#define MOVEMENT_IN_DOMAIN_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <functional>

#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"
#include "./superdropwithgbxindex.hpp"
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

inline void set_gridboxes_superdropletspan(std::vector<GridBox> &gridboxes,
                                           std::vector<SuperdropWithGbxindex> &SDsInGBxs)
{
  for (auto &gbx : gridboxes)
  {
    gbx.set_span(SDsInGBxs);
    
    //gbx.iscorrect_span_for_gbxindex(mdlmaps);
  }
}

inline void exchange_superdroplets_between_gridboxes(std::vector<SuperdropWithGbxindex> &SDsInGBxs,
                                                     std::vector<GridBox> &gridboxes)
/* move superdroplets between gridboxes by changing their associated
gridboxindex if necessary, then (re)sorting SDsInGBxs vector and
updating spans4SDsInGbx for each gridbox */
{
  sort_superdrops_via_gridboxindex(SDsInGBxs);
  set_gridboxes_superdropletspan(gridboxes, SDsInGBxs);
}

template <typename BackwardIdxFunc, typename ForwardIdxFunc>
unsigned int changeindex_ifcoord_outofbounds(const unsigned int index,
                                      const std::pair<double, double> bounds,
                                      const double coord,
                                      const BackwardIdxFunc backwardsidx,
                                      const ForwardIdxFunc forwardsidx)
/* Given bounds = {lowerbound, upperbound} of a gridbox with
index 'gbxindex', function determines if coord is within bounds
of that gridbox. (Note: lower bound exclusive, upper bound inclusive).
If coord not within bounds backwardsidx or forwardsidx function, 
as appropriate, is used to return a neighbouring gridbox's index.
If coord lies within bounds, gbxindex is returned */
{
  if (coord < bounds.first) // lowerbound
  {
    return backwardsidx(index);
  }
  else if (coord >= bounds.second) // upperbound
  {
    return forwardsidx(index);
  }
  else
  {
    return index; // no change to index if coord within bounds
  }
}

#endif // MOVEMENT_IN_DOMAIN_HPP