/*
 * ----- CLEO -----
 * File: microphysicalprocess.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 16th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Microphysical Process Concept as well as
 * helpers for creating structures that obey the
 * concept to model microphysics in SDM, 
 * eg. condensation or collision-coalescence
 * (see ConstTstepProcess struct)
 */

#ifndef MICROPHYSICALPROCESS_HPP
#define MICROPHYSICALPROCESS_HPP

#include <concepts>

template <typename P>
concept MicrophysicalProcess = requires(P p, const unsigned int t)
/* concept for Microphysical Process is all types that
meet requirements (constraints) of these two timstepping
functions ()"on_step" and "next_step") as well as the
constraints on the "run_step" function */
{
  {
    p.next_step(t)
  } -> std::convertible_to<unsigned int>;
  {
    p.on_step(t)
  } -> std::same_as<bool>;
  {
    p.run_step(t)
  } -> std::same_as<void>; 
};

#endif // MICROPHYSICALPROCESS_HPP