/*
 * ----- CLEO -----
 * File: motion.hpp
 * Project: superdrops
 * Created Date: Monday 16th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 17th October 2023
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

template <typename M>
concept Motion = requires(M m, const unsigned int t)
/* concept for superdrop motion is all types that
meet requirements (constraints) of these two timstepping
functions ("on_step" and "next_step") as well as the
constraints on the "update_superdrop_coords" function */
{
  {
    m.next_step(t)
  } -> std::convertible_to<unsigned int>;
  {
    m.on_step(t)
  } -> std::same_as<bool>;
  {
    m.update_superdrop_coords(t)
  } -> std::same_as<void>; 
};

#endif // MOTION_HPP