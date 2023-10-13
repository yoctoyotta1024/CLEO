/*
 * ----- CLEO -----
 * File: gridbox.hpp
 * Project: sdmdomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 13th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Functions and structures related to the CLEO gridboxes
 */


#ifndef GRIDBOX_HPP 
#define GRIDBOX_HPP 

#include <iterator>
#include <algorithm>
#include <span>
#include <string>
#include <stdexcept>
#include <utility>

#include <Kokkos_Core.hpp>

#include "./detectors.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

struct Gridbox
/* each gridbox has unique identifier and contains a
reference to superdroplets in gridbox, alongside the
Gridbox's State (e.g. thermodynamic variables
used for SDM) and detectors for tracking chosen variables */
{
  const unsigned int gbxindex;      // index (unique identifier) of gridbox
  Detectors detectors;              // detectors of various quantities
  std::span<Superdrop> supersingbx; // superdrops in gridbox
  State state;                      // dynamical state of gridbox (e.g. thermodynamics)

  KOKKOS_INLINE_FUNCTION GridBox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~GridBox() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  GridBox(const unsigned int igbxindex)
      : gbxindex(igbxindex) {}
}

struct Gridboxes
{
  // kokkos array (dualview)
};

#endif // GRIDBOX_HPP