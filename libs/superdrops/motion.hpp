/*
 * ----- CLEO -----
 * File: motion.hpp
 * Project: superdrops
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 7th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * concept for motion of superdroplets
 * by changing their spatial coordinates
 */


#ifndef MOTION_HPP
#define MOTION_HPP

#include <concepts>

#include "../cleoconstants.hpp"
#include "./state.hpp"
#include "./superdrop.hpp"

template <typename M, typename GbxMaps>
concept Motion = requires(M m,
                          const unsigned int u,
                          const GbxMaps &gbxmaps,
                          const State &state,
                          Superdrop &sd)
/* concept for superdrop motion is all types that
meet requirements (constraints) of these two timstepping
functions ("on_step" and "next_step") as well as the
constraints on the "update_superdrop_coords" function */
{
  {
    m.next_step(u)
  } -> std::convertible_to<unsigned int>;
  {
    m.on_step(u)
  } -> std::same_as<bool>;
  {
    m.update_superdrop_coords(u, gbxmaps, state, sd)
  } -> std::same_as<void>;
};

struct NullMotion
{
  unsigned int next_step(const unsigned int t_mdl) const
  {
    return LIMITVALUES::uintmax;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return false;
  }

  void update_superdrop_coords(const unsigned int t_mdl) const {}
};
#endif // MOTION_HPP