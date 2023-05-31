// Author: Clara Bayley
// File: "observers.hpp"
/* Observer Concept and structures for
various ways of observing a gridbox of the
superdroplet model, gridbox contained the
thermodynamic state and the vector of
superdroplets' state. Observation is
for example printing some thermodynamics data
to terminal or writing them to a datafile */

#ifndef OBSERVERS_HPP
#define OBSERVERS_HPP

#include <concepts>
#include <iostream>
#include <ios>
#include <iomanip>
#include <stdexcept>
#include <limits>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "../claras_SDconstants.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "superdrop_solver/thermostate.hpp"

namespace dlc = dimless_constants;

template <typename Obs>
concept Observer = requires(Obs obs, const int t, const size_t n,
                            const Kokkos::View<GridBox *> h_gbxs)
/* concept Observer is all types that have an (lvalue) integer
called 'interval' and a function called observe_state() which
take a view of gridboxes as an argument and returns a void type */
{
  {
    obs.interval
  } -> std::same_as<int&>;
  {
    obs.on_step(t)
  } -> std::convertible_to<bool>;
  {
    obs.observe_state(n, h_gbxs)
  } -> std::same_as<void>;
};

template <Observer O1, Observer O2>
struct CombinedObserver
/* combination of two observers is observer 1 followed by observer 2 */
{
  O1 observer1;
  O2 observer2;
  int interval;

  CombinedObserver(const O1 observer1, const O2 observer2)
      : observer1(observer1), observer2(observer2),
        interval(observer1.interval)
  {
    if ((interval != observer1.interval) ||
        (interval != observer2.interval))
    {
      throw std::invalid_argument("observer intervals must be equal");
    }
  }

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gbxs) const
  {
    observer1.observe_state(ngbxs, h_gbxs);
    observer2.observe_state(ngbxs, h_gbxs);
  }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }
};

auto operator>>(const Observer auto o1, const Observer auto o2)
/* define ">>" operator that combines two observers */
{
  return CombinedObserver(o1, o2);
}

struct NullObserver
/* NullObserver does nothing at all
(is defined for a Monoid Structure) */
{
  int interval = std::numeric_limits<int>::max();

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const {}

  bool on_step(const int t) const
  {
    return false;
  }
};

struct PrintObserver
/* this observer prints some details about the
thermodynamic state and superdroplets to terminal */
{
  int interval;            // interval (integer timestep) between observations 
  const int printprec = 4; // precision to print data with

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const;
  /* print time, thermodynamic data (p, temp, qv, qc)
  and total number of superdrops to terminal */

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }
};

#endif // OBSERVERS_HPP