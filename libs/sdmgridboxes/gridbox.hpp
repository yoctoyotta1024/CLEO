// Author: Clara Bayley
// File: "gridbox.hpp"
/* Header file for functions and
structures related to the gridbox
struct */

#ifndef GRIDBOX_HPP
#define GRIDBOX_HPP

#include <vector>
#include <iterator>
#include <algorithm>
#include <span>
#include <string>
#include <stdexcept>

#include "./maps4gridboxes.hpp"
#include "./superdrops_with_gridboxes.hpp"
#include "initialisation/config.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"

struct GridBox
/* gridbox contains vector of superdroplets in grid box,
thermodynamic state temp, pressure, etc. used for SDM,
and index for finding associated grridbox in
coupled thermodynamics */
{
  unsigned int gbxindex; // index / unique identifier of gridbox
  std::span<SuperdropWithGridbox> span4SDsinGBx;
  ThermoState state;

  GridBox(const unsigned int ii,
          const Maps4GridBoxes &mdlmaps,
          std::vector<SuperdropWithGridbox> &SDsInGBxs);

  void set_span(std::vector<SuperdropWithGridbox> &SDsInGBxs);
  /* assumes SDsInGBxs is ordered based on sd_gbxindex
  from lowest to highest. Finds first and last SDWithGBx that has 
  sd_gbxindex matching gbxindex in order to set span4SDsinGBx. */

  void set_statevolume(const std::map<unsigned int, double> &idx2vol);
  /* set dimensionless value for gridbox state's 
  volume using Map4GridBoxes idx2vol map.
  True volume = state.volume * COORD0^3 [m^3] */

  void iscorrect_span_for_gbxindex(const Maps4GridBoxes &mdlmaps);
  /* throw error if the coordinates of the superdroplets
  in the span do not lie within the gridboux boundaries
  given my the gbxindex */

  void iscoord_within_bounds(const std::map<unsigned int,
                                            std::pair<double,
                                                      double>> &idx2bounds,
                             const double coord);
};

std::vector<GridBox> create_gridboxes(const Maps4GridBoxes &mdlmaps,
                                      std::vector<SuperdropWithGridbox> &SDsInGBxs);
/* create domain as a vector of grid boxes such that each grid box
is initialised with a labels from mdlmaps.gbxidxs, and a span of the
superdroplet 'SDsInGbxs', and an (uninitialised) thermodynamic state. */

inline void set_gridboxes_superdropletspan(std::vector<GridBox> &gridboxes,
                                           std::vector<SuperdropWithGridbox> &SDsInGBxs)
{
  for (auto &gbx : gridboxes)
  {
    gbx.set_span(SDsInGBxs);
  }
}

#endif // GRIDBOX_HPP