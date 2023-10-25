/*
 * ----- CLEO -----
 * File: observers.hpp
 * Project: observers
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
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

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

template <typename Obs>
concept Observer = requires(Obs obs, unsigned int t,
                            const viewh_constgbx h_gbxs)
/* concept Observer is all types that have an operator that
has signature of observing functions (see Observer concept) */
{
  {
    obs.before_timestepping(h_gbxs)
  } -> std::same_as<void>;
  {
  obs.next_obs(t)
  } -> std::convertible_to<unsigned int>;
  {
    obs.on_step(t)
  } -> std::same_as<bool>;
  {
    obs.at_start_step(t, h_gbxs)
  } -> std::same_as<void>;
};

template <Observer Obs1, Observer Obs2>
struct CombinedObserver
/* new observer formed from combination
of two Obervers 'a' and 'b' */
{
private:
  Obs1 a;
  Obs2 b;

public:
  CombinedObserver(const Obs1 obs1, const Obs2 obs2)
      : a(obs1), b(obs2) {}

  void before_timestepping(const viewh_constgbx h_gbxs) const
  /* for combination of 2 observers, each
  observer is run sequentially */
  {
    a.before_timestepping(h_gbxs);
    b.before_timestepping(h_gbxs);
  }

  unsigned int next_obs(const unsigned int t_mdl) const
   /* for combination of 2 observers, the next obs
   time is smaller out of the two possible */
  {
    const unsigned int t_a(a.next_obs(t_mdl));
    const unsigned int t_b(b.next_obs(t_mdl));

    return !(t_a < t_b) ? t_b : t_a; // return smaller of two unsigned ints (see std::min)
  }

  bool on_step(const unsigned int t_mdl) const
  /* for combination of 2 observers, a tstep is
  on_step = true when either observer is on_step */
  {
    return a.on_step(t_mdl) || b.on_step(t_mdl);
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  /* for combination of 2 observers, each
  observer is run sequentially */
  {
    a.at_start_step(t_mdl, h_gbxs);
    b.at_start_step(t_mdl, h_gbxs);
  }
};

auto operator>>(const Observer auto obs1,
                const Observer auto obs2)
/* define ">>" operator that combines two observers */
{
  return CombinedObserver(obs1, obs2);
}


struct NullObserver
/* NullObserver does nothing at all
(is defined for a Monoid Structure) */
{
  void before_timestepping(const viewh_constgbx h_gbxs) const {}

  unsigned int next_obs(const unsigned int t_mdl) const
  {
    return LIMITVALUES::uintmax;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return false;
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const {}
};

template <typename O>
concept ObsFuncs = requires(O o, unsigned int t,
                           const viewh_constgbx h_gbxs)
/* concept for all types that can be called used
by ConstTstepObserver for 'do_obs' (in order
to make an Observer type out of a ConstTstepObserver) */
{
  {
    o.before_timestepping(h_gbxs)
  } -> std::same_as<void>;
  {
    o.at_start_step(t, h_gbxs)
  } -> std::same_as<void>;
};

template <ObsFuncs O>
struct ConstTstepObserver
/* this structure is a type that satisfies the concept of
an observer in SDM and has a constant tstep
'interval'. It can be used to create an observer
with a constant timestep and other functionality 
determined by the ObsFunc type 'O' */
{
private:
  unsigned int interval;
  O do_obs;

public:
  ConstTstepObserver(const unsigned int interval, const O o)
      : interval(interval), do_obs(o) {}

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    do_obs.before_timestepping(h_gbxs);
  }

  unsigned int next_obs(const unsigned int t_mdl) const
  {
    return ((t_mdl / interval) + 1) * interval;
  }

  bool on_step(const unsigned int t_mdl) const
  {
    return t_mdl % interval == 0;
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  {
    if (on_step(t_mdl))
    {
      do_obs.at_start_step(t_mdl, h_gbxs);
    }
  }
};

#endif // OBSERVERS_HPP