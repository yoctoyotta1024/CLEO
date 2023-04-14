// Author: Clara Bayley
// File: "intostore_observers.hpp"
/* structures obeying the Observer
concept for various ways of observing
a gridbox which ends up writing data
into a (zarr) store on disk */

#ifndef INTOSTORE_OBSERVERS_HPP
#define INTOSTORE_OBSERVERS_HPP

#include <vector>
#include <string>
#include <stdexcept>

#include "./observers.hpp"
#include "./thermostatestorage.hpp"
#include "./sdattributes_intostore.hpp"
#include "./contigraggedsdstorage.hpp"
#include "./singlevarstorage.hpp"
#include "./sdmomentsstorage.hpp"
#include "sdmgridboxes/gridbox.hpp"

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

class SDMassNthMomentObserver
{
  private:
    const int nth_moment;
    TwoDStorage<double> &zarr;
  
  public:
    SDMassNthMomentObserver(TwoDStorage<double> &zarr,
                            const int nth_moment)
        : nth_moment(nth_moment),
          zarr(zarr)
    {
      if (zarr.get_name() != "massmoment"+std::to_string(nth_moment))
      {
        const std::string errmsg = "name of storage meant for nth "
                                   "moment of SD mass distirbution "
                                   "is not called 'massmoment[n]'";
        throw std::invalid_argument(errmsg);
      }
    }
    
    void observe_state(const std::vector<GridBox> &gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordStorage instance */
    {
      const double n = nth_moment;
      for (auto &gbx : gridboxes)
      {
        double moment = massnthmoment(gbx.span4SDsinGBx, n);
        zarr.value_to_storage(moment);
      }

      ++zarr.nobs;
    }
};

template <typename ContiguousRaggedSDStorage>
class SDsAttributeObserver
{
private:
  ContiguousRaggedSDStorage &zarr;

public:
  SDsAttributeObserver(ContiguousRaggedSDStorage &zarr) : zarr(zarr) {}
  
  void observe_state(const std::vector<GridBox> &gridboxes) const
  /* observe superdroplets by writing their data to contigious
  ragged represented arrays as determined by the 
  ContiguousRaggedSDStorage instance */
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

#endif // INTOSTORE_OBSERVERS_HPP