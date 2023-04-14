// Author: Clara Bayley
// File: "movesuperdropsindomain.hpp"
/* Header file for functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#ifndef MOVEMENT_IN_DOMAIN_HPP
#define MOVEMENT_IN_DOMAIN_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <concepts>

#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"
#include "./superdropwithgbxindex.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/sdmotion.hpp"

/* ----- function called internally ----- */
unsigned int zdown(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int zup(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int xbehind(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int xinfront(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int yleft(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int yright(const Maps4GridBoxes &gbxmaps, const unsigned int index);
/* -------------------------------------- */

template <SdMotion MoveSuperdrop>
class MoveSuperdropsInDomain
{
private:
  const MoveSuperdrop movesd;

  unsigned int update_superdrop_gbxindex(const Maps4GridBoxes &gbxmaps,
                                         const unsigned int gbxindex,
                                         const std::pair<double, double> zbounds,
                                         const std::pair<double, double> xbounds,
                                         const std::pair<double, double> ybounds,
                                         const Superdrop &superdrop) const
  /* For each direction (z, then x, then y), gbxmaps's forward and backward
  get_neighbour functions are passed into changeindex_ifcoord_outofbounds
  along with superdroplet's coord and the gridbox bounds for that direction.
  If coord not within bounds, changeindex_ifcoord_outofbounds is used to
  return a new value of sd_gbxindex via calling the appropriate get_neighbour
  function. After algorithm for z, then x, then y directions are complete,
  resultant sd_gbxindex is returned. */
  {
    unsigned int sd_gbxindex(gbxindex);
    sd_gbxindex = changeindex_ifcoord_outofbounds(gbxmaps,
                                                  zdown, zup,
                                                  zbounds,
                                                  superdrop.coord3,
                                                  sd_gbxindex);

    sd_gbxindex = changeindex_ifcoord_outofbounds(gbxmaps,
                                                  xbehind, xinfront,
                                                  xbounds,
                                                  superdrop.coord1,
                                                  sd_gbxindex);

    sd_gbxindex = changeindex_ifcoord_outofbounds(gbxmaps,
                                                  yleft, yright,
                                                  ybounds,
                                                  superdrop.coord2,
                                                  sd_gbxindex);

    return sd_gbxindex;
}

  void move_superdroplets_between_gridboxes(std::vector<SuperdropWithGbxindex> &SDsInGBxs,
                                            std::vector<GridBox> &gridboxes) const
  /* move superdroplets between gridboxes by changing their associated
  gridboxindex if necessary, then (re)sorting SDsInGBxs vector and
  updating spans4SDsInGbx for each gridbox */
  {
    sort_superdrops_via_gridboxindex(SDsInGBxs);
    set_gridboxes_superdropletspan(gridboxes, SDsInGBxs);
  }

  void set_gridboxes_superdropletspan(std::vector<GridBox> &gridboxes,
                                      std::vector<SuperdropWithGbxindex> &SDsInGBxs) const
  {
    for (auto &gbx : gridboxes)
    {
      gbx.set_span(SDsInGBxs);
      // gbx.iscorrect_span_for_gbxindex(gbxmaps);
    }
  }

  template <typename BackwardIdxFunc, typename ForwardIdxFunc>
  unsigned int changeindex_ifcoord_outofbounds(const Maps4GridBoxes &gbxmaps,
                                               const BackwardIdxFunc backwardsidx,
                                               const ForwardIdxFunc forwardsidx,
                                               const std::pair<double, double> bounds,
                                               const double coord,
                                               const unsigned int sd_gbxindex) const
  /* Given bounds = {lowerbound, upperbound} of a gridbox with
  index 'gbxindex', function determines if coord is within bounds
  of that gridbox. (Note: lower bound inclusive, upper bound exclusive).
  If coord not within bounds backwardsidx or forwardsidx function,
  as appropriate, is used to return a neighbouring gridbox's index.
  If coord lies within bounds, gbxindex is returned. If index is
  already out of domain (ie. value is the maximum unsigned int),
  return out of domain index */
  {
    if (sd_gbxindex == (unsigned int)-1)
    {
      return sd_gbxindex; // sd_gbxindex is out of domain
    }

    if (coord < bounds.first) // lowerbound
    {
      return backwardsidx(gbxmaps, sd_gbxindex);
    }
    else if (coord >= bounds.second) // upperbound
    {
      return forwardsidx(gbxmaps, sd_gbxindex);
    }
    else
    {
      return sd_gbxindex; // no change to index if coord within bounds
    }
  }

public:
  MoveSuperdropsInDomain(const MoveSuperdrop movesd)
      : movesd(movesd) {}

  void move_superdrops_in_domain(const Maps4GridBoxes &gbxmaps,
                                 std::vector<SuperdropWithGbxindex> &SDsInGBxs,
                                 std::vector<GridBox> &gridboxes) const
  /* Move superdroplets that are in gridboxes including exchange
  between gridboxes if necessary. First update superdroplet positions
  according to their motion and then move superdroplets between
  gridboxes by changing their associated gridboxindex as appropriate.
  Final step is (re)sorting SDsInGBxs vector and updating
  spans4SDsInGbx for each gridbox */
  {
    for (auto &gbx : gridboxes)
    {
      const auto zbounds(gbxmaps.get_bounds_z(gbx.gbxindex));
      const auto xbounds(gbxmaps.get_bounds_x(gbx.gbxindex));
      const auto ybounds(gbxmaps.get_bounds_y(gbx.gbxindex));

      for (auto &SDinGBx : gbx.span4SDsinGBx)
      {
        movesd.change_superdroplet_coords(gbx.state, SDinGBx.superdrop);

        SDinGBx.sd_gbxindex = update_superdrop_gbxindex(gbxmaps,
                                                        gbx.gbxindex,
                                                        zbounds, xbounds, ybounds,
                                                        SDinGBx.superdrop);
      }
    }

    move_superdroplets_between_gridboxes(SDsInGBxs, gridboxes);
  }
};

#endif // MOVEMENT_IN_DOMAIN_HPP