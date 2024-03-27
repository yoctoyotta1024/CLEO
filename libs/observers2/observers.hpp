/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: observers.hpp
 * Project: observers2
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 27th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer Concept and related structures for various ways
 * of observing (outputing data from) CLEO.
 * An example of an observer is printing some data
 * from a gridbox's thermostate to the terminal
 */

#ifndef LIBS_OBSERVERS2_OBSERVERS_HPP_
#define LIBS_OBSERVERS2_OBSERVERS_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

/**
 * @brief Concept Observer is all types that have functions for timestepping and
 * at_start_step as constrained here.
 *
 * @tparam Obs The type that satisfies the Observer concept.
 */
template <typename Obs>
concept Observer = requires(Obs obs, unsigned int t, const viewh_constgbx h_gbxs,
                            const viewd_constsupers totsupers, const Gridbox &gbx) {
  { obs.before_timestepping(h_gbxs) } -> std::same_as<void>;
  { obs.after_timestepping() } -> std::same_as<void>;
  { obs.next_obs(t) } -> std::convertible_to<unsigned int>;
  { obs.on_step(t) } -> std::same_as<bool>;
  { obs.at_start_step(t, h_gbxs, totsupers) } -> std::same_as<void>;
  { obs.at_start_step(t, gbx) } -> std::same_as<void>;
};

/**
 * @brief Structure CombinedObserver represents a new observer formed from combination of two
 * Observers 'a' and 'b'.
 *
 * @tparam Obs1 Type satisfying the Observer concept.
 * @tparam Obs2 Type satisfying the Observer concept.
 */
template <Observer Obs1, Observer Obs2>
struct CombinedObserver {
 private:
  Obs1 a; /**< First Observer. */
  Obs2 b; /**< Second Observer. */

 public:
  /**
   * @brief Construct a new CombinedObserver object.
   *
   * @param obs1 First Observer.
   * @param obs2 Second Observer.
   */
  CombinedObserver(const Obs1 obs1, const Obs2 obs2) : a(obs1), b(obs2) {}

  /* for combination of 2 observers, each
  observer is run sequentially */
  void before_timestepping(const viewh_constgbx h_gbxs) const {
    a.before_timestepping(h_gbxs);
    b.before_timestepping(h_gbxs);
  }

  /* for combination of 2 observers, each
  observer is run sequentially */
  void after_timestepping() const {
    a.after_timestepping();
    b.after_timestepping();
  }

  /* for combination of 2 observers, the next obs
  time is smaller out of the two possible */
  unsigned int next_obs(const unsigned int t_mdl) const {
    const auto t_a = a.next_obs(t_mdl);
    const auto t_b = b.next_obs(t_mdl);

    return !(t_a < t_b) ? t_b : t_a;  // return smaller of two unsigned ints (see std::min)
  }

  /* for combination of 2 observers, a tstep is
  on_step = true when either observer is on_step */
  bool on_step(const unsigned int t_mdl) const { return a.on_step(t_mdl) || b.on_step(t_mdl); }

  /* for combination of 2 observers, each
  observer is run sequentially */
  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {
    a.at_start_step(t_mdl, h_gbxs, totsupers);
    b.at_start_step(t_mdl, h_gbxs, totsupers);
  }

  /* for combination of 2 observers, each
  observer is run sequentially */
  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {
    a.at_start_step(t_mdl, gbx);
    b.at_start_step(t_mdl, gbx);
  }
};

/**
 * @brief Overloaded operator >> to combine two Observers.
 *
 * @param obs1 First Observer.
 * @param obs2 Second Observer.
 * @return CombinedObserver<Obs1, Obs2> Combined Observer.
 */
auto operator>>(const Observer auto obs1, const Observer auto obs2) {
  return CombinedObserver(obs1, obs2);
}

/**
 * @brief Structure NullObserver does nothing at all.
 *
 * NullObserver defined for completion of Observer's Monoid Set.
 *
 */
struct NullObserver {
  void before_timestepping(const viewh_constgbx h_gbxs) const {}

  void after_timestepping() const {}

  unsigned int next_obs(const unsigned int t_mdl) const { return LIMITVALUES::uintmax; }

  bool on_step(const unsigned int t_mdl) const { return false; }

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {}

  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {}
};

/**
 * @brief Concept ObsFuncs for all types that can be called used by ConstTstepObserver for
 * observation functions.
 *
 * Type in ConstTstepObserver obeying ObsFunc makes it possible for ConstTstepObserver to obey
 * Observer concept.
 *
 * @tparam O Type that satisfies the ObsFuncs concept.
 */
template <typename O>
concept ObsFuncs = requires(O o, unsigned int t, const viewh_constgbx h_gbxs,
                            const viewd_constsupers totsupers, const Gridbox &gbx) {
  { o.before_timestepping(h_gbxs) } -> std::same_as<void>;
  { o.after_timestepping() } -> std::same_as<void>;
  { o.at_start_step(t, h_gbxs, totsupers) } -> std::same_as<void>;
  { o.at_start_step(t, gbx) } -> std::same_as<void>;
};

/**
 * @brief Structure ConstTstepObserver represents a type that satisfies the concept of an observer
 * with a constant timestep interval between observations.
 *
 * Struct can be used to create an observer with a constant timestep and with observation
 * functionality as determined by the ObsFunc type 'O'
 *
 * @tparam O Type that satisfies the ObsFuncs concept.
 */
template <ObsFuncs O>
struct ConstTstepObserver {
 private:
  unsigned int interval; /**< interval between observations. */
  O do_obs;              /**< Observation functionality. */

 public:
  /**
   * @brief Construct a new ConstTstepObserver object.
   *
   * @param interval Timestep interval.
   * @param o Observer.
   */
  ConstTstepObserver(const unsigned int interval, const O o) : interval(interval), do_obs(o) {}

  void before_timestepping(const viewh_constgbx h_gbxs) const {
    do_obs.before_timestepping(h_gbxs);
  }

  void after_timestepping() const { do_obs.after_timestepping(); }

  unsigned int next_obs(const unsigned int t_mdl) const {
    return ((t_mdl / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_mdl) const { return t_mdl % interval == 0; }

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {
    if (on_step(t_mdl)) {
      do_obs.at_start_step(t_mdl, h_gbxs, totsupers);
    }
  }

  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {
    if (on_step(t_mdl)) {
      do_obs.at_start_step(t_mdl, gbx);
    }
  }
};

#endif  // LIBS_OBSERVERS2_OBSERVERS_HPP_
