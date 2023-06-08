// Author: Clara Bayley
// File: "movesuperdropsindomain.hpp"
/* Header file for functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#ifndef MOVESUPERDROPSINDOMAIN_HPP 
#define MOVESUPERDROPSINDOMAIN_HPP 

#include <map>
#include <utility>
#include <stdexcept>
#include <concepts>

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>

#include "./gridbox.hpp"
#include "./maps4gridboxes.hpp"
#include "./cartesianneighbours.hpp"
#include "./superdropwithgbxindex.hpp"
#include "./sdmotion.hpp"
#include "./detectors.hpp"
#include "superdrop_solver/superdrop.hpp"

/* ----- function called internally -----
  to update superdrop coord and return
  sd_gbxindex of neighbouring gridbox in
  one of 6 particular directions */
unsigned int zdown(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop);
unsigned int zup(const Maps4GridBoxes &gbxmaps,
                 const unsigned int index,
                 Superdrop &superdrop);
unsigned int xbehind(const Maps4GridBoxes &gbxmaps,
                     const unsigned int index,
                     Superdrop &superdrop);
unsigned int xinfront(const Maps4GridBoxes &gbxmaps,
                      const unsigned int index,
                      Superdrop &superdrop);
unsigned int yleft(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop);
unsigned int yright(const Maps4GridBoxes &gbxmaps,
                    const unsigned int index,
                    Superdrop &superdrop);
/* -------------------------------------- */

template <SdMotion MoveSuperdrop>
class MoveSuperdropsInDomain
{
private:
  const MoveSuperdrop movesd;

  void move_superdrops_in_domain(const Maps4GridBoxes &gbxmaps,
                                 Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs,
                                 Kokkos::vector<GridBox> &gridboxes) const
  /* Move superdroplets in gridboxes used movesd and then move them
  between gridboxes if necessary. First update superdroplet positions
  according to their motion and then move superdroplets between
  gridboxes by changing their associated gridboxindex as appropriate.
  Final step is (re)sorting SDsInGBxs vector and updating
  spans4SDsInGbx for each gridbox */
  {
    for (auto &gbx : gridboxes)
    {
      const auto ii(gbx.gbxindex);
      const auto zbds(gbxmaps.get_bounds_z(ii));
      const auto xbds(gbxmaps.get_bounds_x(ii));
      const auto ybds(gbxmaps.get_bounds_y(ii));
      const auto area(gbxmaps.get_area(ii));

      for (auto &SDinGBx : gbx.span4SDsinGBx)
      {
        auto &drop(SDinGBx.superdrop);
        movesd.change_superdroplet_coords(gbxmaps, gbx, drop);

        gbx.detectors -> detect_precipitation(area, drop);

        SDinGBx.sd_gbxindex = update_superdrop_gbxindex(gbxmaps, ii,
                                                        zbds, xbds, ybds,
                                                        drop);
      }
    }

    move_superdroplets_between_gridboxes(gbxmaps, SDsInGBxs, gridboxes);
  }

  unsigned int update_superdrop_gbxindex(const Maps4GridBoxes &gbxmaps,
                                         const unsigned int gbxindex,
                                         std::pair<double, double> zbounds,
                                         std::pair<double, double> xbounds,
                                         std::pair<double, double> ybounds,
                                         Superdrop &superdrop) const
  /* For each direction (z, then x, then y), gbxmaps's forward and backward
  get_neighbour functions are passed into update_superdrop_ifneighbour
  along with superdroplet and the gridbox bounds for that direction.
  (If coord not within bounds, update_superdrop_ifneighbour should call
  appropriate get_neighbour function to update the superdroplet's 
  sd_gbxindex (and possibly other attributes if desired). After algorithm
  for z, then x, then y directions are complete, resultant sd_gbxindex
  is returned. */
  {
    unsigned int current_gbxindex(gbxindex);

    const auto get_zbounds = [&](const auto ii)
    { return gbxmaps.get_bounds_z(ii); };
    const auto get_coord3 = [&]()
    { return superdrop.coord3; };
    current_gbxindex = update_ifneighbour(gbxmaps, zdown, zup,
                                          get_zbounds, get_coord3,
                                          current_gbxindex, zbounds,
                                          superdrop);

    const auto get_xbounds = [&](const unsigned int ii)
    { return gbxmaps.get_bounds_x(ii);};
    const auto get_coord1 = [&]()
    { return superdrop.coord1; };
    current_gbxindex = update_ifneighbour(gbxmaps, xbehind, xinfront,
                                          get_xbounds, get_coord1,
                                          current_gbxindex, xbounds,
                                          superdrop);

    const auto get_ybounds = [&](const unsigned int ii)
    { return gbxmaps.get_bounds_y(ii);};
    const auto get_coord2 = [&]()
    { return superdrop.coord2; };
    current_gbxindex = update_ifneighbour(gbxmaps, yleft, yright,
                                          get_ybounds, get_coord2,
                                          current_gbxindex, ybounds,
                                          superdrop);

    return current_gbxindex;
  }

  template <typename BackwardIdxFunc, typename ForwardIdxFunc,
            typename GetBounds, typename GetSdCoord>
  unsigned int update_ifneighbour(const Maps4GridBoxes &gbxmaps,
                                  const BackwardIdxFunc backwards_neighbour,
                                  const ForwardIdxFunc forwards_neighbour,
                                  const GetBounds get_bounds,
                                  const GetSdCoord get_sdcoord,
                                  unsigned int current_gbxindex,
                                  std::pair<double, double> bounds,
                                  Superdrop &superdrop) const
  /* For a given direction, pass {lower, upper} bounds into
  update_superdrop_ifneighbour to get updated gbxindex and superdrop
  (e.g. if superdroplet's coord fromsdcoord function lies outside of
  bounds given gbxbounds using current_gbxindex). Repeat until
  superdroplet coord is within the bounds given by the current_gbxindex,
  or until superdrop leaves domain. */
  {
    double coord(get_sdcoord());

    /* loop while coord is within domain but not within bounds
    break if coord is out of domain (or within bounds). */
    while ((current_gbxindex != OUTOFDOMAIN()) &&
           ((coord < bounds.first) | (coord >= bounds.second)))
    {
      current_gbxindex = update_superdrop_ifneighbour(gbxmaps,
                                                      backwards_neighbour,
                                                      forwards_neighbour,
                                                      current_gbxindex,
                                                      bounds, coord,
                                                      superdrop);
      bounds = get_bounds(current_gbxindex);
      coord = get_sdcoord();
    }

    return current_gbxindex;
  }

  template <typename BackwardIdxFunc, typename ForwardIdxFunc>
  unsigned int update_superdrop_ifneighbour(const Maps4GridBoxes &gbxmaps,
                                            const BackwardIdxFunc backwards_neighbour,
                                            const ForwardIdxFunc forwards_neighbour,
                                            const unsigned int current_gbxindex,
                                            const std::pair<double, double> bounds,
                                            const double coord,
                                            Superdrop &superdrop) const
  /* Given bounds = {lowerbound, upperbound} of a gridbox with
  index 'gbxindex', function determines if coord is within bounds
  of that gridbox. (Note: lower bound inclusive, upper bound exclusive).
  If coord not within bounds backwardsidx or forwardsidx function,
  as appropriate, is used to return a neighbouring gridbox's index.
  If coord lies within bounds, gbxindex is returned. If index is
  already out of domain (ie. value is the maximum unsigned int),
  return out of domain index */
  {
    if (current_gbxindex == OUTOFDOMAIN())
    {
      return current_gbxindex; // ignore SDs whose sd_gbxindex is already out of domain
    }

    if (coord < bounds.first) // lowerbound
    {
      return backwards_neighbour(gbxmaps, current_gbxindex, superdrop);
    }
    else if (coord >= bounds.second) // upperbound
    {
      return forwards_neighbour(gbxmaps, current_gbxindex, superdrop);
    }
    else
    {
      return current_gbxindex; // no change to index if coord within bounds
    }
  }

  void move_superdroplets_between_gridboxes(const Maps4GridBoxes &gbxmaps,
                                            Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs,
                                            Kokkos::vector<GridBox> &gridboxes) const
  /* move superdroplets between gridboxes by changing their associated
  gridboxindex if necessary, then (re)sorting SDsInGBxs vector and
  updating spans4SDsInGbx for each gridbox */
  {
    sort_superdrops_via_gridboxindex(SDsInGBxs);
    set_gridboxes_superdropletspan(gbxmaps, gridboxes, SDsInGBxs);
  }

  void set_gridboxes_superdropletspan(const Maps4GridBoxes &gbxmaps,
                                      Kokkos::vector<GridBox> &gridboxes,
                                      Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs) const
  {
    for (auto &gbx : gridboxes)
    {
      gbx.set_span(SDsInGBxs);
      // gbx.iscorrect_span_for_gbxindex(gbxmaps); // (expensive!) optional test to raise error if SDspan isn't consistent with gbxindex 
    }
  }

public:
  MoveSuperdropsInDomain(const MoveSuperdrop movesd)
      : movesd(movesd) {}

  int next_step(const int currenttimestep) const
  {
    return movesd.next_move(currenttimestep);
  }

  void run_step(const int currenttimestep,
           const Maps4GridBoxes &gbxmaps,
           Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs,
           Kokkos::vector<GridBox> &gridboxes) const
  {
    if (movesd.on_move(currenttimestep))
    {
      move_superdrops_in_domain(gbxmaps, SDsInGBxs, gridboxes);
    }
  }
};

#endif // MOVESUPERDROPSINDOMAIN_HPP 