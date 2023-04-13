
// Author: Clara Bayley
// File: "movesuperdropsindomain.cpp"
/* Implementation of functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#include "./movesuperdropsindomain.hpp"

/* ----- function called internally ----- */
unsigned int zdown(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int zup(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int xbehind(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int xinfront(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int yleft(const Maps4GridBoxes &gbxmaps, const unsigned int index);
unsigned int yright(const Maps4GridBoxes &gbxmaps, const unsigned int index);
/* -------------------------------------- */

unsigned int MoveSuperdropsInDomain::update_superdrop_gbxindex(
    const Maps4GridBoxes &gbxmaps,
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

unsigned int zdown(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index)
{
  return gbxmaps.get_neighbour_zdown(index);
};

unsigned int zup(const Maps4GridBoxes &gbxmaps,
                 const unsigned int index)
{
  return gbxmaps.get_neighbour_zup(index);
};

unsigned int xbehind(const Maps4GridBoxes &gbxmaps,
                     const unsigned int index)
{
  return gbxmaps.get_neighbour_xbehind(index);
};

unsigned int xinfront(const Maps4GridBoxes &gbxmaps,
                      const unsigned int index)
{
  return gbxmaps.get_neighbour_xinfront(index);
};

unsigned int yleft(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index)
{
  return gbxmaps.get_neighbour_yleft(index);
};

unsigned int yright(const Maps4GridBoxes &gbxmaps,
                    const unsigned int index)
{
  return gbxmaps.get_neighbour_yright(index);
};