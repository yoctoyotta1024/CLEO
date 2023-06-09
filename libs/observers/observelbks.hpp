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

#include "../claras_SDconstants.hpp"
#include "sdmgridboxes/logbooks.hpp"
#include "sdmgridboxes/gridbox.hpp"

namespace dlc = dimless_constants;

template <typename OL>
concept ObserveLbks = requires(OL o, const DetectorLogbooks lbks)
/* concept ObserveLbks is all types that have an operator that
has signature of observe_logbooks() function (see Observer concept)
ie. which takes a logbooks struct as argument and returns void */
{
  {
  o(lbks)
  } -> std::same_as<void>;
  {
    o.prepare(lbks)
  } -> std::same_as<void>;
};

template <ObserveLbks Ol1, ObserveLbks Ol2>
class CombinedObserveLbks
/* combination of two types obeying the ObserveLbks
concept is 'ol1' followed by 'ol2' (resultant
combination also obeys ObserveLbks concept) */
{
private:
  Ol1 o1;
  Ol2 o2;

public:
  CombinedObserveLbks(const Ol1 o1, const Ol2 o2)
      : o1(o1), o2(o2) {}

  void operator()(const DetectorLogbooks &lbks) const
  {
    o1(lbks);
    o2(lbks);
  }

  void prepare(const DetectorLogbooks &lbks) const
  {
    o1.prepare(lbks);
    o2.prepare(lbks);
  }
};

auto operator>>(const ObserveLbks auto o1,
                const ObserveLbks auto o2)
/* define ">>" operator that combines
two observe logbooks types */
{
  return CombinedObserveLbks(o1, o2);
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

  ObsLbks obslbks;

public:
  ConstIntervalLbksObserver(const int interval,
                            const ObsLbks observe_logbooks)
      : interval(interval), obslbks(observe_logbooks) {}

  int get_interval() const { return interval; }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  void prepare(const DetectorLogbooks &lbks) const
  {
    obslbks.prepare(lbks);
  }

  void observe_logbooks(const DetectorLogbooks &lbks) const
  {
    obslbks(lbks);
  }

  void observe_gridboxes(const size_t ngbxs,
                         const Kokkos::View<GridBox *> h_gridboxes) const {}

  void observe(const size_t ngbxs,
               const Kokkos::View<GridBox *> h_gridboxes,
               const DetectorLogbooks &lbks) const
  {
    observe_logbooks(ngbxs, h_gridboxes);
  }
};


struct PrintLogbooks
/* satisfies ObserveLbks concept and
prints out details about logbooks */
{
  void printprecip(
      const std::shared_ptr<Logbook<double>> &logbook) const
  {
    constexpr int printprec(6); // precision to print data with
    
    std::cout << "dimless surfpp = ";

    for (size_t idx = 0; idx < logbook  -> get_size(); ++idx)
    {
      std::cout << std::scientific
                << std::setprecision(printprec)
                << logbook->get_entry(idx) << ", ";
    }
    std::cout << '\n';
  }

  void prepare(const DetectorLogbooks &lbks) const {}

  void operator()(const DetectorLogbooks &logbooks) const
  {
    printprecip(logbooks.surfpp); 
  }
};

#endif // OBSERVELBKS_HPP