/*
 * ----- CLEO -----
 * File: predcorrmotion.hpp
 * Project: gridboxes
 * Created Date: Tuesday 19th December 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 19th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * -----
 * File Description:
 * Generic struct satisfying Motion concept for 
 * a superdroplet using predictor-corrector
 * method to update a superdroplet's coordinates and
 * updating gbx according to templated functions
 */


#ifndef PREDCORRMOTION_HPP 
#define PREDCORRMOTION_HPP 

#include <functional>
#include <cassert>

#include <Kokkos_Core.hpp>

#include "./predcorr.hpp"
#include "superdrops/superdrop.hpp"
#include "superdrops/terminalvelocity.hpp"

template <GridboxMaps GbxMaps>
struct ChangeIfCoord3Nghbr
{
  CartesianDomainBoundaryCoord3 if_domainboundary_coord3;

  KOKKOS_INLINE_FUNCTION unsigned int
  change_to_backwards_coord3nghbr(const unsigned int idx,
                                  const GbxMaps &gbxmaps,
                                  Superdrop &drop) const
  /* function to update superdrop and return gbxindex of
  neighbouring gridbox in backwards coord3 direction
  (including case when superdrop exceeds domain boundary) */
  {
    auto nghbr = (unsigned int)gbxmaps.coord3backward(idx);

    nghbr = if_domainboundary_coord3.backwards(idx, nghbr, gbxmaps, drop);

    drop.set_sdgbxindex(nghbr);

    return nghbr; // gbxindex of coord3 backwards neighbour
  };

  KOKKOS_INLINE_FUNCTION unsigned int
  change_to_forwards_coord3nghbr(const unsigned int idx,
                                 const GbxMaps &gbxmaps,
                                 Superdrop &drop) const
  /* function to update superdrop and return gbxindex of
  neighbouring gridbox in forwards coord3 direction
  (including case when superdrop exceeds domain boundary) */
  {
    auto nghbr = (unsigned int)gbxmaps.coord3forward(idx);

    nghbr = if_domainboundary_coord3.forwards(idx, nghbr, gbxmaps, drop);

    drop.set_sdgbxindex(nghbr);

    return nghbr; // gbxindex of coord3 forwards neighbour
  };

  KOKKOS_INLINE_FUNCTION unsigned int
  operator()(const int flag,
             const GbxMaps &gbxmaps,
             unsigned int idx,
             Superdrop &drop) const
  /* change superdrop if it should move to neighbouring
  gridbox in coord3 direction. Function alters superdrop
  and returns new value of gbxindex if flag != 0:
  if flag = 1 idx updated to backwards neighbour gbxindex.
  if flag = 2 idx updated to forwards neighbour gbxindex.
  Note: backwards/forwards functions may change the
  superdroplet's attributes e.g. if it leaves the domain. */
  {
    switch (flag)
    {
    case 1:
      idx = change_to_backwards_coord3nghbr(idx, gbxmaps, drop);
      break;
    case 2:
      idx = change_to_forwards_coord3nghbr(idx, gbxmaps, drop);
      break;
    }
    return idx;
  }
};

template <GridboxMaps GbxMaps, typename Flag>
struct ChangeIfNeighbour
/* wrapper of functions to update superdrop sdgbxindex
and possibly other attributes if superdrop should
move to neighbouring gridbox e.g. see usage in
PredCorrMotion. Struct has three functions, one for each
direction (coord3 = z, coord1 = x, coord2 = y). For each,
the superdrop's coord is compared to gridbox bounds given by gbxmaps
for the current gbxindex 'idx'. If superdrop coord lies outside
bounds, forward or backward neighbour functions are called to
update sdgbxindex (and possibly other superdrop attributes) */
{
private:
  Flag flag_sdgbxindex;
  ChangeIfCoord3Nghbr<GbxMaps> change_if_coord3nghbr;
  ChangeIfCoord1Nghbr<GbxMaps> change_if_coord1nghbr;
  ChangeIfCoord2Nghbr<GbxMaps> change_if_coord2nghbr;

public:
  ChangeIfNeighbour(const Flag i_flag_sdgbxindex,
                    ChangeIfCoord3Nghbr<GbxMaps> i_change_if_coord3nghbr)
      : flag_sdgbxindex(i_flag_sdgbxindex),
        change_if_coord3nghbr(i_change_if_coord3nghbr) {}

  KOKKOS_INLINE_FUNCTION unsigned int
  coord3(const GbxMaps &gbxmaps,
         unsigned int idx,
         Superdrop &drop) const
  /* change superdrop if it should move to neighbouring
  gridbox in coord3 direction. Function alters superdrop
  and returns new value of gbxindex based on flag */
  {
    const auto flag = flag_sdgbxindex(idx,
                                      gbxmaps.coord3bounds(idx),
                                      drop.get_coord3()); // if value != 0 idx needs to change

    idx = change_if_coord3nghbr(flag, gbxmaps, idx, drop);

    return idx;
  }

  KOKKOS_INLINE_FUNCTION unsigned int
  coord1(const CartesianMaps &gbxmaps,
         unsigned int idx,
         Superdrop &drop) const
  /* change superdrop if it should move to neighbouring
  gridbox in coord1 direction. Function alters superdrop
  and returns new value of gbxindex based on flag */
  {
    const auto flag = flag_sdgbxindex(idx,
                                      gbxmaps.coord1bounds(idx),
                                      drop.get_coord1()); // if value != 0 idx needs to change

    idx = change_if_coord1nghbr(flag, gbxmaps, idx, drop);

    return idx;
  }

  KOKKOS_INLINE_FUNCTION unsigned int
  coord2(const CartesianMaps &gbxmaps,
         unsigned int idx,
         Superdrop &drop) const
  /* change superdrop if it should move to neighbouring
  gridbox in coord1 direction. Function alters superdrop
  and returns new value of gbxindex based on flag */
  {
    const auto flag = flag_sdgbxindex(idx,
                                      gbxmaps.coord2bounds(idx),
                                      drop.get_coord2()); // if value != 0 idx needs to change

    idx = change_if_coord2nghbr(flag, gbxmaps, idx, drop);
    
    return idx;
  }
};

template <GridboxMaps GbxMaps, VelocityFormula TV,
          typename CheckBounds, typename ChangeIfNghbr>
struct PredCorrMotion
/* satisfies motion concept for motion of a superdroplet
using a predictor-corrector method to update a superdroplet's
coordinates and then updating it's sdgbxindex using
the appropriate templated type */
{
  const unsigned int interval; // integer timestep for movement
  PredCorr<GbxMaps, TV> superdrop_coords;
  CheckBounds check_bounds;
  ChangeIfNghbr change_if_nghbr;

  PredCorrMotion(const unsigned int motionstep,
                 const std::function<double(unsigned int)> int2time,
                 const TV i_terminalv,
                 CheckBounds i_check_bounds,
                 ChangeIfNghbr i_change_if_nghbr)
      : interval(motionstep),
        superdrop_coords(interval, int2time, i_terminalv),
        check_bounds(i_check_bounds),
        change_if_nghbr(i_change_if_nghbr) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_sdm) const
  {
    return ((t_sdm / interval) + 1) * interval;
  }

  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int t_sdm) const
  {
    return t_sdm % interval == 0;
  }

  KOKKOS_INLINE_FUNCTION void
  superdrop_gbx(const unsigned int gbxindex,
                const CartesianMaps &gbxmaps,
                Superdrop &drop) const
  /* function satisfies requirements of
  "superdrop_gbx" in the motion concept to update a
  superdroplet if it should move between gridboxes.
  For each direction (coord3, then coord1, then coord2),
  superdrop and idx may be changed if superdrop coord
  lies outside bounds of gridbox in that direction */
  {
    unsigned int idx = change_if_nghbr.coord3(gbxmaps, gbxindex, drop);
    check_bounds(idx, gbxmaps.coord3bounds(idx), drop.get_coord3());

    idx = change_if_nghbr.coord1(gbxmaps, idx, drop);
    check_bounds(idx, gbxmaps.coord1bounds(idx), drop.get_coord1());

    idx = change_if_nghbr.coord2(gbxmaps, idx, drop);
    check_bounds(idx, gbxmaps.coord2bounds(idx), drop.get_coord2());

    assert((drop.get_sdgbxindex() == idx) &&
           "sdgbxindex not concordant with supposed idx");
  }
};

/* -----  ----- TODO: move functions below to .cpp file ----- ----- */

#endif // PREDCORRMOTION_HPP  