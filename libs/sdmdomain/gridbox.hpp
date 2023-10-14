/*
 * ----- CLEO -----
 * File: gridbox.hpp
 * Project: sdmdomain
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 14th October 2023
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

struct SupersInGridbox
/* Reference to a chunk of memory
(e.g. through std::span of Kokkos::subview)
containing super-droplets in a Gridbox */
{
  size_t num;               // number of superdrops in gridbox
  std::span<Superdrop> sds; // reference to superdrops in gridbox

  KOKKOS_INLINE_FUNCTION SupersInGridbox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~SupersInGridbox() = default; // Kokkos requirement for a (dual)View
};

struct Gridbox
/* each gridbox has unique identifier and contains a
reference to superdroplets in gridbox, alongside the
Gridbox's State (e.g. thermodynamic variables
used for SDM) and detectors for tracking chosen variables */
{
private:
  unsigned int gbxindex;             // index (unique identifier) of gridbox
public:
  Detectors detectors;               // detectors of various quantities
  SupersInGridbox supersingbx;       // reference to superdrops associated with gridbox
  State state;                       // dynamical state of gridbox (e.g. thermodynamics)

  KOKKOS_INLINE_FUNCTION Gridbox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Gridbox() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  Gridbox(const unsigned int igbxindex,
          const double ivolume)
      : gbxindex(igbxindex),
        detectors(),
        supersingbx(),
        state(ivolume) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int get_gbxindex() {return gbxindex;}
};

#endif // GRIDBOX_HPP