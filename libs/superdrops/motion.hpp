/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: motion.hpp
 * Project: superdrops
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * concept for motion of superdroplets by changing their spatial coordinates
 */

#ifndef LIBS_SUPERDROPS_MOTION_HPP_
#define LIBS_SUPERDROPS_MOTION_HPP_

#include <concepts>

#include "../cleoconstants.hpp"
#include "state.hpp"
#include "superdrop.hpp"

/* concept for superdrop motion is all types that
meet requirements (constraints) of these two timstepping
functions ("on_step" and "next_step") as well as the
constraints on the "superdrop_coords" function */
template <typename M, typename GbxMaps>
concept Motion =
    requires(M m, const unsigned int u, const GbxMaps &gbxmaps, const State &state, Superdrop &sd) {
      { m.next_step(u) } -> std::convertible_to<unsigned int>;
      { m.on_step(u) } -> std::same_as<bool>;
      { m.superdrop_coords(u, gbxmaps, state, sd) } -> std::same_as<void>;
      { m.superdrop_gbx(u, gbxmaps, sd) } -> std::same_as<void>;
    };

struct NullMotion {
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int t_mdl) const { return LIMITVALUES::uintmax; }

  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int t_mdl) const { return false; }

  template <typename GbxMaps>
  KOKKOS_INLINE_FUNCTION void superdrop_coords(const unsigned int gbxindex, const GbxMaps &gbxmaps,
                                               const State &state, Superdrop &drop) const {}

  template <typename GbxMaps>
  KOKKOS_INLINE_FUNCTION void superdrop_gbx(const unsigned int gbxindex, const GbxMaps &gbxmaps,
                                            Superdrop &drop) const {}
};

#endif  // LIBS_SUPERDROPS_MOTION_HPP_
