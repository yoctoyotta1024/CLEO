
// Author: Clara Bayley
// File: "movesuperdropsindomain.cpp"
/* Implementation of functions related to
moving superdroplets (both updating their
coords and moving them between gridboxes) */

#include "./movesuperdropsindomain.hpp"

unsigned int zdown_periodic(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_zdown(index);
};


unsigned int zdown_finite(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_zdown(index);
};

unsigned int zup(const Maps4GridBoxes &gbxmaps,
                 const unsigned int index,
                 Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_zup(index);
};

unsigned int xbehind(const Maps4GridBoxes &gbxmaps,
                     const unsigned int index,
                     Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_xbehind(index);
};

unsigned int xinfront(const Maps4GridBoxes &gbxmaps,
                      const unsigned int index,
                      Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_xinfront(index);
};

unsigned int yleft(const Maps4GridBoxes &gbxmaps,
                   const unsigned int index,
                   Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_yleft(index);
};

unsigned int yright(const Maps4GridBoxes &gbxmaps,
                    const unsigned int index,
                    Superdrop &superdrop)
{
  return gbxmaps.get_neighbour_yright(index);
};