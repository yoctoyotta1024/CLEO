/*
 * ----- CLEO -----
 * File: observers.hpp
 * Project: observers
 * Created Date: Friday 13th October 2023
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
 * Observer Concept and related structures for various ways
 * of observing (outputing data from) CLEO.
 * An example of an observer is printing some data
 * from a gridbox's thermostate to the terminal
 */


#ifndef OBSERVERS_HPP
#define OBSERVERS_HPP

#include <concepts>

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "sdmdomain/gridbox.hpp"

template <typename Obs>
concept Observer = requires(Obs obs, unsigned int t,
                            viewh_constgbx h_gbxs)
/* concept Observer is all types that have an operator that
has signature of observing functions (see Observer concept) */
{
  {
    obs.get_obsstep()
  } -> std::convertible_to<unsigned int>;
  {
    obs.on_step(t)
  } -> std::same_as<bool>;
  {
    obs.observe_startstep(t, h_gbxs)
  } -> std::same_as<void>;
};

#endif // OBSERVERS_HPP