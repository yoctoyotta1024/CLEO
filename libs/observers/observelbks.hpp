// Author: Clara Bayley
// File: "observelbks.hpp"
/* ObserveLbks Concept and related 
structures for various ways of observing 
logbooks of the superdroplet model. An
example of an observe logbook type may be
something that writes a the data from a
logbook to an array in a zarr storage system */

#ifndef OBSERVELBKS_HPP
#define OBSERVELBKS_HPP

#include <concepts>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

template <typename L>
concept ObserveLbks = requires(L l, const std::vector<int> lbks)
/* concept ObserveLbks is all types that have an operator that
has signature of observe_logbooks() function (see Observer concept)
ie. which takes a logbooks struct as argument and returns void */
{
  {
  obs(lbks)
  } -> std::same_as<void>;
};

template <ObserveLbks Ol1, ObserveLbks Ol2>
class CombinedObserveLbks
/* combination of two ObserveGridBox types
is 'og1' followed by 'og2' */
{
private:
  Ol o1;
  Ol o2;

public:
  CombinedObserveLbks(const Ol o1, const Ol o2)
      : o1(o1), o2(o2) {}

  void operator()(const std::vector<int> lbks) const
  {
    o1(lbks);
    o2(lbks);
  }
};

auto operator>>(const ObserveLbks auto o1,
                const ObserveLbks auto o2)
/* define ">>" operator that combines
two observe logbooks types */
{
  return CombinedObserveGBxs(o1, o2);
}

template <ObserveLbks ObsLbks>
class ConstIntervalLbksObserver
/* struct satifying the Observer concept
that has constant time-step 'interval'
between obseration of gridboxes and
takes no action during observe_logbooks */
{
private:
  const int interval; // interval (integer timestep) between observations

  ObsLbks observe_logbooks;

public:
  ConstIntervalLbksObserver(const int interval, const ObsLbks obslbks)
      : interval(interval), observe_logbooks(obslbks) {}

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gridboxes) const {}

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes,
               const std::vector<int> lbks) const
  {
    observe_logbooks(ngbxs, h_gridboxes);
  }
};

#endif // OBSERVELBKS_HPP