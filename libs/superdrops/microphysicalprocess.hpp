/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: microphysicalprocess.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 14th December 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Microphysical Process Concept as well as
 * helpers for creating structures that obey the
 * concept to model microphysics in SDM,
 * eg. condensation or collision-coalescence
 * (see ConstTstepProcess struct)
 */

#ifndef LIBS_SUPERDROPS_MICROPHYSICALPROCESS_HPP_
#define LIBS_SUPERDROPS_MICROPHYSICALPROCESS_HPP_

#include <concepts>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include "../cleoconstants.hpp"
#include "./kokkosaliases_sd.hpp"
#include "./state.hpp"
#include "./superdrop.hpp"

/* concept for Microphysical Process is all types that
meet requirements (constraints) of these two timstepping
functions ()"on_step" and "next_step") as well as the
constraints on the "run_step" function */
template <typename P>
concept MicrophysicalProcess = requires(P p, const TeamMember &tm, const unsigned int t,
                                        subviewd_supers supers, State &state, GenRandomPool gp) {
  { p.next_step(t) } -> std::convertible_to<unsigned int>;
  { p.on_step(t) } -> std::same_as<bool>;
  { p.run_step(tm, t, supers, state, gp) } -> std::convertible_to<subviewd_supers>;
};

template <MicrophysicalProcess Microphys1, MicrophysicalProcess Microphys2>
struct CombinedMicrophysicalProcess {
 private:
  Microphys1 a;
  Microphys2 b;

 public:
  CombinedMicrophysicalProcess(const Microphys1 a, const Microphys2 b) : a(a), b(b) {}

  /* for combination of 2 microphysical processes,
  the next timestep is smaller out of the two possible */
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int subt) const {
    const auto t_a = a.next_step(subt);
    const auto t_b = b.next_step(subt);

    return !(t_a < t_b) ? t_b : t_a;  // return smaller of two unsigned ints (see std::min)
  }

  /* for combination of 2 microphysical processes,
  a tstep is on_step = true if either individual
  process is on_step = true */
  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int subt) const { return a.on_step(subt) || b.on_step(subt); }

  /* for combination of 2 proceses, each process
  is called sequentially */
  KOKKOS_INLINE_FUNCTION subviewd_supers run_step(const TeamMember &team_member,
                                                  const unsigned int subt, subviewd_supers supers,
                                                  State &state, GenRandomPool genpool) const {
    supers = a.run_step(team_member, subt, supers, state, genpool);
    supers = b.run_step(team_member, subt, supers, state, genpool);
    return supers;
  }
};

/* define ">>" operator that combines
two Superdroplet Model Microphysical Processes */
auto operator>>(const MicrophysicalProcess auto a, const MicrophysicalProcess auto b) {
  return CombinedMicrophysicalProcess(a, b);
}

/* NullProcess does nothing at all
(is defined for a Monoid Structure) */
struct NullMicrophysicalProcess {
  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int subt) const { return LIMITVALUES::uintmax; }

  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int subt) const { return false; }

  KOKKOS_INLINE_FUNCTION subviewd_supers run_step(const TeamMember &team_member,
                                                  const unsigned int subt, subviewd_supers supers,
                                                  State &state, GenRandomPool genpool) const {
    return supers;
  }
};

/* concept for all (function-like) types
(ie. types that can be called with some arguments)
that can be called by the run_step function in
ConstTstepMicrophysics (see below) */
template <typename F>
concept MicrophysicsFunc = requires(F f, const TeamMember &tm, const unsigned int subt,
                                    subviewd_supers supers, State &state, GenRandomPool gp) {
  { f(tm, subt, supers, state, gp) } -> std::convertible_to<subviewd_supers>;
};

/* this structure is a type that satisfies the concept of
microphysical process in SDM and has a constant tstep
'interval'. It can be used to create a microphysical
processes with a constant timestep and microphysics
determined by the MicrophysicsFunc type 'F' */
template <MicrophysicsFunc F>
struct ConstTstepMicrophysics {
 private:
  unsigned int interval;
  F do_microphysics;

 public:
  ConstTstepMicrophysics(const unsigned int interval, const F f)
      : interval(interval), do_microphysics(f) {}

  KOKKOS_INLINE_FUNCTION
  unsigned int next_step(const unsigned int subt) const {
    return ((subt / interval) + 1) * interval;
  }

  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int subt) const { return subt % interval == 0; }

  KOKKOS_INLINE_FUNCTION subviewd_supers run_step(const TeamMember &team_member,
                                                  const unsigned int subt, subviewd_supers supers,
                                                  State &state, GenRandomPool genpool) const {
    if (on_step(subt)) {
      supers = do_microphysics(team_member, subt, supers, state, genpool);
    }

    return supers;
  }
};

#endif  // LIBS_SUPERDROPS_MICROPHYSICALPROCESS_HPP_
