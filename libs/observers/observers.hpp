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
#include <vector>
#include <iostream>
#include <ios>
#include <iomanip>

#include "sdmgridboxes/gridbox.hpp"
#include "superdrop_solver/thermostate.hpp"

template <typename Obs>
concept Observer = requires(Obs obs,
                            const std::vector<GridBox> &gboxes)
/* concept Observer is all types that have a function
called observe_state() which take a gridbox type as
argument and returns a void type */
{
  {
    obs.observe_state(gboxes)
    } -> std::same_as<void>;
};

template <Observer O1, Observer O2>
struct CombinedObserver
/* combination of two observers is observer 1 followed by observer 2 */
{
  O1 observer1;
  O2 observer2;

  CombinedObserver(const O1 observer1, const O2 observer2)
      : observer1(observer1), observer2(observer2) {}

  void observe_state(const std::vector<GridBox> &gboxes) const
  {
    observer1.observe_state(gboxes);
    observer2.observe_state(gboxes);
  }
};

auto operator>>(const Observer auto o1, const Observer auto o2)
/* define ">>" operator that combines two observers */
{
  return CombinedObserver{o1, o2};
}

struct NullObserver
/* NullObserver does nothing at all
(is defined for a Monoid Structure) */
{
  void observe_state(const std::vector<GridBox> &gridboxes) const {}
};

struct PrintObserver
/* this observer prints some details about the
thermodynamic state and superdroplets to terminal */
{
  const int printprec = 4; // precision to print data with

  void observe_state(const std::vector<GridBox> &gridboxes) const;
  /* print time, thermodynamic data (p, temp, qv, qc)
  and total number of superdrops to terminal */
};

#endif // OBSERVERS_HPP