// Author: Clara Bayley
// File: "gridbox.hpp"
/* Header file for functions and
structures related to the gridbox
struct */

#ifndef GRIDBOX_HPP
#define GRIDBOX_HPP

#include <iterator>
#include <algorithm>
#include <span>
#include <string>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>

#include "../claras_SDconstants.hpp"
#include "./maps4gridboxes.hpp"
#include "./superdropwithgbxindex.hpp"
#include "./logbooks.hpp"
#include "./detectors.hpp"
#include "initialisation/config.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"
#include "superdrop_solver/thermodynamic_equations.hpp"

namespace dlc = dimless_constants;

struct GridBox
/* gridbox contains vector of superdroplets in grid box,
thermodynamic state temp, pressure, etc. used for SDM,
and index for finding associated grridbox in
coupled thermodynamics */
{
  unsigned int gbxindex; // index (unique identifier) of gridbox
  std::span<SuperdropWithGbxindex> span4SDsinGBx;
  ThermoState state;
  std::unique_ptr<Detectors> detectors;

  KOKKOS_INLINE_FUNCTION GridBox() = default; // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~GridBox() = default; // Kokkos requirement for a (dual)View

  KOKKOS_FUNCTION
  GridBox(const unsigned int ii,
          const Maps4GridBoxes &gbxmaps,
          const InstallDetectors &dtrs,
          Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs);
  /* Volume in Thermostate set using Map4GridBoxes
  idx2vol map (via get_volume function). Other ThermoState variables
  are default behaviour initialised. */

  KOKKOS_FUNCTION
  void set_span(Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs);
  /* assumes SDsInGBxs is ordered based on sd_gbxindex
  from lowest to highest. Finds first and last SDWithGBx that has 
  sd_gbxindex matching gbxindex in order to set span4SDsinGBx. */

  KOKKOS_FUNCTION
  void print_statevolume();
  /* print's dimensionless value for gridbox state's 
  volume. Also prints true volume = volume * COORD0^3 [m^3] */

  KOKKOS_FUNCTION
  void iscorrect_span_for_gbxindex(const Maps4GridBoxes &gbxmaps);
  /* throw error if the coordinates of the superdroplets
  in the span do not lie within the gridboux boundaries
  given my the gbxindex */

  KOKKOS_FUNCTION
  void iscoord_within_bounds(const std::pair<double, double> bounds,
                             const double coord);
};

KOKKOS_FUNCTION Kokkos::vector<GridBox>
create_gridboxes(const Maps4GridBoxes &gbxmaps,
                 Kokkos::vector<SuperdropWithGbxindex> &SDsInGBxs);
/* create domain as a vector of grid boxes such that each grid box
is initialised with a labels from gbxmaps.gbxidxs, and a span of the
superdroplet 'SDsInGbxs', and an (uninitialised) thermodynamic state. */

KOKKOS_FUNCTION
void set_superdroplets_to_wetradius(Kokkos::vector<GridBox> &gridboxes);
/* for each gridbox, set the radius of each superdroplet (SD) to
whichever is larger out of their dry radius or equlibrium wet radius
(given the relative humidity (s_ratio) and temperature of the gridbox).
If relh > maxrelh = 0.95, set each SD's radius to their
equilibrium radius at relh = maxrelh = 0.95 */

#endif // GRIDBOX_HPP