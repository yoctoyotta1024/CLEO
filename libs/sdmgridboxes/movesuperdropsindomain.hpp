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
      const auto zbounds(gbxmaps.get_bounds_z(gbx.gbxindex));
      const auto xbounds(gbxmaps.get_bounds_x(gbx.gbxindex));
      const auto ybounds(gbxmaps.get_bounds_y(gbx.gbxindex));

      for (auto &SDinGBx : gbx.span4SDsinGBx)
      {
        movesd.change_superdroplet_coords(gbxmaps, gbx, SDinGBx.superdrop);

        gbx.detector.precip(zbounds, SDinGBx.superdrop);

        SDinGBx.sd_gbxindex = update_superdrop_gbxindex(gbxmaps,
                                                        gbx.gbxindex,
                                                        zbounds, xbounds, ybounds,
                                                        SDinGBx.superdrop);
      }
    }

    move_superdroplets_between_gridboxes(gbxmaps, SDsInGBxs, gridboxes);
  }

  unsigned int update_superdrop_gbxindex(const Maps4GridBoxes &gbxmaps,
                                         const unsigned int gbxindex,
                                         const std::pair<double, double> zbounds,
                                         const std::pair<double, double> xbounds,
                                         const std::pair<double, double> ybounds,
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
    unsigned int sd_gbxindex(gbxindex);
    sd_gbxindex = update_superdrop_ifneighbour(gbxmaps,
                                              zdown, zup,
                                              zbounds,
                                              superdrop.coord3,
                                              sd_gbxindex,
                                              superdrop);

    sd_gbxindex = update_superdrop_ifneighbour(gbxmaps,
                                        xbehind, xinfront,
                                        xbounds,
                                        superdrop.coord1,
                                        sd_gbxindex,
                                        superdrop);

    sd_gbxindex = update_superdrop_ifneighbour(gbxmaps,
                                                yleft, yright,
                                                ybounds,
                                                superdrop.coord2,
                                                sd_gbxindex,
                                                superdrop);

    return sd_gbxindex;
  }

  template <typename BackwardIdxFunc, typename ForwardIdxFunc>
  unsigned int update_superdrop_ifneighbour(const Maps4GridBoxes &gbxmaps,
                                            const BackwardIdxFunc backwards_neighbour,
                                            const ForwardIdxFunc forwards_neighbour,
                                            const std::pair<double, double> bounds,
                                            const double coord,
                                            const unsigned int sd_gbxindex,
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
    if (sd_gbxindex == std::numeric_limits<unsigned int>::max())
    {
      return sd_gbxindex; // ignore SDs whose sd_gbxindex is already out of domain
    }

    if (coord < bounds.first) // lowerbound
    {
      return backwards_neighbour(gbxmaps, sd_gbxindex, superdrop);
    }
    else if (coord >= bounds.second) // upperbound
    {
      return forwards_neighbour(gbxmaps, sd_gbxindex, superdrop);
    }
    else
    {
      return sd_gbxindex; // no change to index if coord within bounds
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