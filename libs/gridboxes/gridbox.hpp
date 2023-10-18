/*
 * ----- CLEO -----
 * File: gridbox.hpp
 * Project: gridboxes
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
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

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

#include "./detectors.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

struct Gridbox
/* each gridbox has unique identifier and contains a
reference to superdroplets in gridbox, alongside the
Gridbox's State (e.g. thermodynamic variables
used for SDM) and detectors for tracking chosen variables */
{
private:
  struct SupersInGbx
  /* References to identify the chunk of memory
  containing super-droplets occupying a given Gridbox
  (e.g. through std::span or Kokkos::subview) */
  {
  private:
    Kokkos::pair<size_t, size_t> pos; // position in view of (first, last) superdrop that occupies gridbox

  public:
    KOKKOS_INLINE_FUNCTION SupersInGbx() = default;  // Kokkos requirement for a (dual)View
    KOKKOS_INLINE_FUNCTION ~SupersInGbx() = default; // Kokkos requirement for a (dual)View

    KOKKOS_INLINE_FUNCTION
    SupersInGbx(const Kokkos::pair<size_t, size_t> ipos)
        : pos(ipos) {}

    template <class view_type>
    Kokkos::Subview<view_type, Kokkos::pair<size_t, size_t>>
    operator()(view_type supers)
    /* returns subview from view of superdrops
    refering to superdrops which occupy given gridbox */
    {
      return Kokkos::subview(supers, pos.first, pos.second);
    }
  };

public:
  struct Gbxindex
  /* struct containing gridbox index
  and its generator struct */
  {
    unsigned int value;
    class Gen
    {
    public:
      Gbxindex next() { return {_id++}; }

    private:
      unsigned int _id = 0;
    };
  };

  Gbxindex gbxindex;    // index (unique identifier) of gridbox
  State state;          // dynamical state of gridbox (e.g. thermodynamics)
  SupersInGbx sdsingbx; // reference(s) to superdrops occupying gridbox
  Detectors detectors;  // detectors of various quantities

  KOKKOS_INLINE_FUNCTION Gridbox() = default;  // Kokkos requirement for a (dual)View
  KOKKOS_INLINE_FUNCTION ~Gridbox() = default; // Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  Gridbox(const Gbxindex igbxindex,
          const State istate)
      : gbxindex(igbxindex),
        state(istate),
        sdsingbx(Kokkos::make_pair<size_t, size_t>(0,0)),
        detectors() {}

  KOKKOS_INLINE_FUNCTION
  auto get_gbxindex() {return gbxindex.value;}
};

#endif // GRIDBOX_HPP