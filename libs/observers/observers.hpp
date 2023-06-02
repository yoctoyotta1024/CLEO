// Author: Clara Bayley
// File: "observers.hpp"
/* Observer Concept and related structures 
for various ways of observing gridboxes 
and logbooks of the superdroplet model. 
An example of an observer is printing some data
from a gridbox's thermostate to the terminal */

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

#include "./observegbxs.hpp"
#include "./observelbks.hpp"
#include "../claras_SDconstants.hpp"
#include "sdmgridboxes/gridbox.hpp"
#include "sdmgridboxes/logbooks.hpp"
#include "superdrop_solver/thermostate.hpp"

namespace dlc = dimless_constants;

template <typename Obs>
concept Observer = requires(Obs obs, const int t, const size_t n,
                            const Kokkos::View<GridBox *> h_gbxs,
                            const DetectorLogbooks lbks)
/* concept Observer is all types that have
functions with signatures and return types
required for observing gridboxes and logbooks */
{
  {
    obs.observe(n, h_gbxs, lbks)
  } -> std::same_as<void>;
  {
    obs.observe_gridboxes(n, h_gbxs)
  } -> std::same_as<void>;
  {
    obs.observe_logbooks(lbks)
  } -> std::same_as<void>;
  {
    obs.get_interval()
  } -> std::convertible_to<int>;
  {
    obs.on_step(t)
  } -> std::convertible_to<bool>;
  {
    obs.prepare()
  } -> std::same_as<void>;
};

template <Observer Obs1, Observer Obs2>
class CombinedObserver
/* combination of two observers 'obs1'
and 'obs2' as long as 'obs1' and 'obs2'
have return same value from get_interval() */
{
private:
  const int interval; // interval (integer timestep) between observations
  Obs1 o1;
  Obs2 o2;

public:
  CombinedObserver(const Obs1 observer1, const Obs2 observer2)
      : interval(observer1.get_interval()),
        o1(observer1), o2(observer2)
  {
    const auto intvl = get_interval();
    if ((intvl != o1.get_interval()) || (intvl != o2.get_interval()))
    {
      throw std::invalid_argument("observer intervals must be equal");
    }
  }

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }
  
  void prepare() const
  {
    o1.prepare();
    o2.prepare();
  }

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gbxs) const
  {
    o1.observe_gridboxes(ngbxs, h_gbxs);
    o2.observe_gridboxes(ngbxs, h_gbxs);
  }

  void observe_logbooks(const DetectorLogbooks &lbks) const
  {
    o1.observe_logbooks(lbks);
    o2.observe_logbooks(lbks);
  }

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gbxs,
               const DetectorLogbooks &lbks) const
  {
    observe_gridboxes(ngbxs, h_gbxs);
    observe_logbooks(lbks);
  }
};

auto operator>>(const Observer auto obs1,
                const Observer auto obs2)
/* define ">>" operator that combines two observers */
{
  return CombinedObserver(obs1, obs2);
}

struct NullObserver
/* NullObserver (also a NullObserveGBxs)
does nothing at all. It is defined for
completion of a Monoid Structure */
{
  void operator()(const size_t ngbxs,
                  const Kokkos::View<GridBox *> h_gbxs) const {}

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gbxs) const {}

  void observe_logbooks(const DetectorLogbooks &lbks) const {}

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gbxs,
               const DetectorLogbooks &lbks) const {}

  void prepare() const {}

  int get_interval() { return std::numeric_limits<int>::max(); }

  bool on_step(const int t) const
  {
    return false;
  }
};

template <ObserveGBxs ObsGBxs, ObserveLbks ObsLbks>
class ConstIntervalObserver
/* struct satifying the Observer concept
that has constant time-step 'interval'
between obserations of gridboxes and
logbooks */
{
private:
  const int interval; // interval (integer timestep) between observations

  ObsGBxs obsgbxs;
  ObsLbks obslbks;

public:
  ConstIntervalObserver(const int interval,
                        const ObsGBxs observe_gridboxes,
                        const ObsLbks observe_logbooks)
      : interval(interval), obsgbxs(observe_gridboxes),
        obslbks(observe_logbooks) {}

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  void prepare() const
  {
    obsgbxs.prepare();
    obslbks.prepare();
  }

  void observe_gridboxes(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes) const
  {
    obsgbxs(ngbxs, h_gridboxes);
  }

   void observe_logbooks(const DetectorLogbooks &logbooks) const
  {
    obslbks(logbooks);
  }

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes,
               const DetectorLogbooks &logbooks) const
  {
    observe_gridboxes(ngbxs, h_gridboxes);
    observe_logbooks(logbooks);
  }
};

struct PrintObserver
/* satisfies Observer concept and
prints out details about logbooks 
and gridboxes' thermodynamic states
and superdroplets */
{
  const int interval;      // interval (integer timestep) between observations

  PrintObserver(const int obsstep) : interval(obsstep) {}

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  void prepare() const {}

  void observe_logbooks(const DetectorLogbooks &logbooks) const;

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gridboxes) const;
  /* print time, thermodynamic data (p, temp, qv, qc)
  and total number of superdrops to terminal */

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes,
               const DetectorLogbooks &logbooks) const
  {
    observe_gridboxes(ngbxs, h_gridboxes);
  }
};

#endif // OBSERVERS_HPP