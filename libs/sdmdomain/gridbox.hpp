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

#include <span>

#include <Kokkos_Core.hpp>

#include "./detectors.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

struct SuperdropsInGridbox
/* Reference to a chunk of memory
(e.g. through std::span of Kokkos::subview)
containing super-droplets in a Gridbox */
{
private:
  Kokkos::pair<size_t, size_t> pos; // position in view of (first, last) superdrop that occupies gridbox
public:
  KOKKOS_INLINE_FUNCTION SuperdropsInGridbox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~SuperdropsInGridbox() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  SuperdropsInGridbox(const Kokkos::pair<size_t, size_t> ipos)
      : pos(ipos) {}

  template <class view_type>
  Kokkos::Subview<view_type, Kokkos::pair<size_t, size_t>>
  operator()(view_type supers)
  /* returns subview from view of superdrops
  for superdrops which occupy given gridbox */
  {
    return Kokkos::subview(supers, pos.first, pos.second);
  }
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
  State state;                       // dynamical state of gridbox (e.g. thermodynamics)
  SuperdropsInGridbox sdsingbx;       // reference to superdrops associated with gridbox
  Detectors detectors;               // detectors of various quantities

  KOKKOS_INLINE_FUNCTION Gridbox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Gridbox() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  Gridbox(const unsigned int igbxindex,
          const double ivolume,
          const Kokkos::pair<size_t, size_t> ipos)
      : gbxindex(igbxindex),
        state(ivolume),
        sdsingbx(ipos),
        detectors() {}

  KOKKOS_INLINE_FUNCTION
  unsigned int get_gbxindex() {return gbxindex;}
};

#endif // GRIDBOX_HPP