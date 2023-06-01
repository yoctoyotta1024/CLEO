// Author: Clara Bayley
// File: "observers.hpp"
/* Observer and ObserveGridBox Concepts
and related structures for various ways of
observing gridboxes and logbooks of the
superdroplet model. An example of an 
observer is printing some data from a
gridbox's thermostate to the terminal */

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
                            const Kokkos::View<GridBox *> h_gbxs,
                            const std::vector<int> lgbks)
/* concept Observer is all types that have
functions with signatures and return types
required for observing gridboxes and logbooks */
{
  {
    obs.observe(n, h_gbxs, lgbks)
  } -> std::same_as<void>;
  {
    obs.observe_gridboxes(n, h_gbxs)
  } -> std::same_as<void>;
  {
    obs.observe_logbooks(lgbks)
  } -> std::same_as<void>;
  {
    obs.get_interval()
  } -> std::convertible_to<int>;
  {
    obs.on_step(t)
  } -> std::convertible_to<bool>;
};

template <typename Og>
concept ObserveGBxs = requires(Og og, const int t, const size_t n,
                               const Kokkos::View<GridBox *> h_gbxs,
                               const std::vector<int> lgbks)
/* concept ObserveGBxs is all types that have an operator that 
has signature of observe_gridboxes() function in Observer concept
ie. which takes a size_t type and a view of gridboxes as an
argument and returns void */
{
  {
    og(n, h_gbxs)
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

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gbxs) const
  {
    o1.observe_gridboxes(ngbxs, h_gbxs);
    o2.observe_gridboxes(ngbxs, h_gbxs);
  }

  void observe_logbooks(const std::vector<int> lgbks) const
  {
    o1.observe_logbooks(lgbks);
    o2.observe_logbooks(lgbks);
  }

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gbxs,
               const std::vector<int> lgbks) const
  {
    observe_gridboxes(ngbxs, h_gbxs);
    observe_logbooks(lgbks);
  }
};

auto operator>>(const Observer auto obs1,
                const Observer auto obs2)
/* define ">>" operator that combines two observers */
{
  return CombinedObserver(obs1, obs2);
}

template <ObserveGBxs Og1, ObserveGBxs Og2>
class CombinedObserveGBxs
/* combination of two ObserveGridBox types
is 'og1' followed by 'og2' */
{
private:
  Og1 og1;
  Og2 og2;

public:
  CombinedObserveGBxs(const Og1 og1, const Og2 og2)
      : og1(og1), og2(og2) {}

  void operator()(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gbxs) const
  {
    og1(ngbxs, h_gbxs);
    og2(ngbxs, h_gbxs);
  }
};

auto operator>>(const ObserveGBxs auto og1,
                const ObserveGBxs auto og2)
/* define ">>" operator that combines two observe gridbox types */
{
  return CombinedObserveGBxs(og1, og2);
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

  void observe_logbooks(const std::vector<int> lgbks) const {}

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gbxs,
               const std::vector<int> lgbks) const {}

  int get_interval() { return std::numeric_limits<int>::max(); }

  bool on_step(const int t) const
  {
    return false;
  }
};

template <ObserveGBxs ObsGBxs>
class ConstIntervalGBxObserver
/* struct satifying the Observer concept
that has constant time-step 'interval'
between obseration of gridboxes and
takes no action during observe_logbooks */
{
private:
  const int interval; // interval (integer timestep) between observations

  ObsGBxs obsgbxs;

public:
  ConstIntervalGBxObserver(const int interval,
                           const ObsGBxs obsgbxs)
      : interval(interval), obsgbxs(obsgbxs) {}

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gridboxes) const
  {
    obsgbxs(ngbxs, h_gridboxes);
  }

  void observe_logbooks(const std::vector<int> lgbks) const {}

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes,
               const std::vector<int> lgbks) const
  {
    obsgbxs(ngbxs, h_gridboxes);
  }
};

struct PrintObserver
/* satisfies Observer concept and
prints out details about gridboxes'
thermodynamic states and superdroplets */
{
  const int interval;      // interval (integer timestep) between observations
  const int printprec = 4; // precision to print data with

  PrintObserver(const int obsstep) : interval(obsstep) {}

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  void observe_logbooks(const std::vector<int> lgbks) const {}

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gridboxes) const;
  /* print time, thermodynamic data (p, temp, qv, qc)
  and total number of superdrops to terminal */

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes,
               const std::vector<int> lgbks) const
  {
    observe_gridboxes(ngbxs, h_gridboxes);
  }
};

#endif // OBSERVERS_HPP