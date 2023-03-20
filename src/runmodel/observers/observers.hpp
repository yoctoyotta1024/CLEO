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
#include <fstream>
#include <string>
#include <span>
#include <iostream>
#include <ios>
#include <iomanip>

#include "./observer_superdropletattributes.hpp"
#include "./observer_thermostate.hpp"
#include "./observer_singlevariable.hpp"
#include "../gridbox.hpp"
#include "superdrop_solver/superdrop.hpp"
#include "superdrop_solver/thermostate.hpp"

template <typename Obs>
concept Observer = requires(Obs obs, const std::vector<GridBox> &gboxes)
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

  CombinedObserver(O1 observer1, O2 observer2)
      : observer1(observer1), observer2(observer2) {}

  void observe_state(const std::vector<GridBox> &gboxes) const
  {
    observer1.observe_state(gboxes);
    observer2.observe_state(gboxes);
  }
};

auto operator>>(Observer auto o1, Observer auto o2)
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

class ThermoStateObserver
{
  private:
    ThermoStateStorage &zarr;
  
  public:
    ThermoStateObserver(ThermoStateStorage &zarr) : zarr(zarr) {}
    
    void observe_state(const std::vector<GridBox> &gridboxes) const
    /* observe thermostate by writing it to arrays
    as determined by the ThermoStateStorage instance */
    {
      for (auto &gbx : gridboxes)
      {
        zarr.thermodata_to_storage(gbx.state);
      }

      ++zarr.nobs;
    }
};

template <typename ContiguousRaggedSuperdropStorage>
class SDsAttributeObserver
{
private:
  ContiguousRaggedSuperdropStorage &zarr;

public:
  SDsAttributeObserver(ContiguousRaggedSuperdropStorage &zarr) : zarr(zarr) {}
  
  void observe_state(const std::vector<GridBox> &gridboxes) const
  /* observe superdroplets by writing their data to contigious
  ragged represented arrays as determined by the 
  ContiguousRaggedSuperdropStorage instance */
  {
    size_t nsupers(0);
    for (auto &gbx : gridboxes)
    {
      for (auto &SDinGBx : gbx.span4SDsinGBx)
      {
        zarr.data_to_contigraggedarray(SDinGBx.superdrop);
        ++nsupers;
      }
    }
    zarr.contigraggedarray_count(nsupers);
  }
};

class TimeObserver
{
  private:
    CoordStorage<double> &zarr;
  
  public:
    TimeObserver(CoordStorage<double> &zarr) : zarr(zarr)
    {
      if (zarr.get_name() != "time")
      {
        const std::string errmsg = "name of storage meant for "
                                   "timeobserver is not called 'time'";
        throw std::invalid_argument(errmsg);
      }
    }
    
    void observe_state(const std::vector<GridBox> &gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordStorage instance */
    {
      auto gbx = gridboxes[0];
      zarr.value_to_storage(gbx.state.time);
    }
};

class GridBoxIndexObserver
{
  private:
    CoordStorage<unsigned int> &zarr;
  
  public:
    GridBoxIndexObserver(CoordStorage<unsigned int> &zarr) : zarr(zarr)
    {
      if (zarr.get_name() != "gbxindex")
      {
        const std::string errmsg = "name of storage meant for gridbox index"
                                   " observer is not called 'gbxindex'";
        throw std::invalid_argument(errmsg);
      }
    }
    
    void observe_state(const std::vector<GridBox> &gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordStorage instance */
    {
      if (zarr.get_ndata() == 0)
      {
        for (auto &gbx : gridboxes)
        {
          zarr.value_to_storage(gbx.gbxindex);
        }
      }
    }
};

class NsupersPerGridBoxObserver
{
  private:
    TwoDStorage<size_t> &zarr;
  
  public:
    NsupersPerGridBoxObserver(TwoDStorage<size_t> &zarr) : zarr(zarr)
    {
      if (zarr.get_name() != "nsupers")
      {
        const std::string errmsg = "name of storage meant for no. superdroplets"
                                   " per gridbox observer is not called 'nsupers'";
        throw std::invalid_argument(errmsg);
      }
    }
    
    void observe_state(const std::vector<GridBox> &gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordStorage instance */
    {
      for (auto &gbx : gridboxes)
      {
        size_t nsupers = gbx.span4SDsinGBx.size();
        zarr.value_to_storage(nsupers);
      }

      ++zarr.nobs;
    }
};

#endif // OBSERVERS_HPP