/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: gridbox.hpp
 * Project: gridboxes
 * Created Date: Wednesday 8th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Functions and structures related to the CLEO gridboxes
 */

#ifndef LIBS_GRIDBOXES_GRIDBOX_HPP_
#define LIBS_GRIDBOXES_GRIDBOX_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "gridboxes/gbxindex.hpp"
#include "gridboxes/supersingbx.hpp"
#include "superdrops/kokkosaliases_sd.hpp"
#include "superdrops/state.hpp"
#include "superdrops/superdrop.hpp"

/* each gridbox has unique identifier and contains a
reference to superdroplets in gridbox, alongside the
Gridbox's State (e.g. thermodynamic variables
used for SDM) */
struct Gridbox {
  Gbxindex gbxindex;        // index (unique identifier) of gridbox
  State state;              // dynamical state of gridbox (e.g. thermodynamics)
  SupersInGbx supersingbx;  // reference(s) to superdrops occupying gridbox

  Gridbox() = default;   // Kokkos requirement for a (dual)View
  ~Gridbox() = default;  // Kokkos requirement for a (dual)View

  /* assumes supers view (or subview) already sorted via sdgbxindex. Constructor works
  outside of parallelism */
  Gridbox(const Gbxindex igbxindex, const State istate, const viewd_supers totsupers)
      : gbxindex(igbxindex), state(istate), supersingbx(totsupers, gbxindex.value) {}

  /* assumes supers view (or subview) already sorted via sdgbxindex. Constructor works within
  parallel team policy on host given member 'team_member' */
  Gridbox(const Gbxindex igbxindex, const State istate, const viewd_supers totsupers,
          const kkpair_size_t irefs)
      : gbxindex(igbxindex), state(istate), supersingbx(totsupers, gbxindex.value, irefs) {}

  KOKKOS_INLINE_FUNCTION
  auto get_gbxindex() const { return gbxindex.value; }

  viewd_constsupers domain_totsupers_readonly() const {
    return supersingbx.domain_totsupers_readonly();
  }

  size_t domain_totnsupers() const { return supersingbx.domain_totnsupers(); }
};

#endif  // LIBS_GRIDBOXES_GRIDBOX_HPP_
