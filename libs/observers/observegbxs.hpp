// Author: Clara Bayley
// File: "observegbxs.hpp"
/* ObserveGBxs Concept and related 
structures for various ways of observing 
gridboxes of the superdroplet model. An
example of an observe gridbox type may be
something that writes a gridboxes' thermostate
to an array in a zarr storage system */

#ifndef OBSERVEGBXS_HPP
#define OBSERVEGBXS_HPP

#include <concepts>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "sdmgridboxes/gridbox.hpp"

template <typename Og>
concept ObserveGBxs = requires(Og og, const size_t n,
                               const Kokkos::View<GridBox *> h_gbxs)
/* concept ObserveGBxs is all types that have an operator that
has signature of observe_gridboxes() function (see Observer concept)
ie. which takes a size_t type and a view of gridboxes as an
argument and returns void */
{
  {
    og(n, h_gbxs)
  } -> std::same_as<void>;
};

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
/* define ">>" operator that combines
two observe gridbox types */
{
  return CombinedObserveGBxs(og1, og2);
}

template <ObserveGBxs ObsGBxs>
class ConstIntervalGBxsObserver
/* struct satifying the Observer concept
that has constant time-step 'interval'
between obseration of gridboxes and
takes no action during observe_logbooks */
{
private:
  const int interval; // interval (integer timestep) between observations

  ObsGBxs observe_gridboxes;

public:
  ConstIntervalGBxsObserver(const int interval, const ObsGBxs obsgbxs)
      : interval(interval), observe_gridboxes(obsgbxs) {}

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  void observe_logbooks(const std::vector<int> lbks) const {}

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes,
               const std::vector<int> lbks) const
  {
    observe_gridboxes(ngbxs, h_gridboxes);
  }
};

#endif // OBSERVEGBXS_HPP