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
called 'interval' and a function called observe_gridboxes() which
take a view of gridboxes as an argument and returns a void type */
{
  {
    obs.observe_gridboxes(n, h_gbxs)
  } -> std::same_as<void>;
  {
    obs.get_interval()
  } -> std::convertible_to<int>;
  {
    obs.on_step(t)
  } -> std::convertible_to<bool>;
};

class ConstIntervalStep
{
private:
  const int interval; // interval (integer timestep) between observations

public:
  ConstIntervalStep(const int interval) : interval(interval) {}

  bool operator()(const int t) const
  /* on_step boolean function */
  {
    return t % interval == 0;
  }

  int get_interval() const { return interval; }
};

template <Observer O1, Observer O2>
class CombinedObserver
/* combination of two observers is observer
'obs1' followed by observer 'obs2' as long as
'obs1' and 'obs2' have same onstep interval */
{
private:
  O1 o1;
  O2 o2;

public:
  ConstIntervalStep on_step;

  CombinedObserver(const O1 observer1, const O2 observer2)
      : o1(observer1), o2(observer2),
        on_step(observer1.get_interval())
  {
    const auto intvl = get_interval();
    if ((intvl != o1.get_interval()) || (intvl != o2.get_interval()))
    {
      throw std::invalid_argument("observer intervals must be equal");
    }
  }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gbxs) const
  {
    o1.observe_gridboxes(ngbxs, h_gbxs);
    o2.observe_gridboxes(ngbxs, h_gbxs);
  }

  int get_interval() const { return on_step.get_interval(); }
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
  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const {}

  int get_interval() { return std::numeric_limits<int>::max(); }

  bool on_step(const int t) const
  {
    return false;
  }
};

struct PrintObserver
/* this observer prints some details about the
thermodynamic state and superdroplets to terminal */
{
  ConstIntervalStep on_step;
  const int printprec = 4; // precision to print data with

  PrintObserver(const int obsstep) : on_step(obsstep) {}

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const;
  /* print time, thermodynamic data (p, temp, qv, qc)
  and total number of superdrops to terminal */

  int get_interval() const { return on_step.get_interval(); }
};

#endif // OBSERVERS_HPP