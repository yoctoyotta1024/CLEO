// Author: Clara Bayley
// File: "intostore_observers.hpp"
/* structures obeying the Observer
concept for various ways of observing
a gridbox which ends up writing data
into a (zarr) store on disk */

#ifndef INTOSTORE_OBSERVERS_HPP
#define INTOSTORE_OBSERVERS_HPP

#include <string>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

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
    
    void observe_state(const Kokkos::View<GridBox*> h_gridboxes) const
    /* observe thermostate by writing it to arrays
    as determined by the ThermoStateStorage instance */
    {
      const size_t Ngrid = h_gridboxes.size();
      for (size_t ii(0); ii<Ngrid; ++ii)
      {
        zarr.thermodata_to_storage(h_gridboxes(ii).state);
      }

      ++zarr.nobs;
    }
};

class TimeObserver
{
  private:
    CoordinateStorage<double> &zarr;
  
  public:
    TimeObserver(CoordinateStorage<double> &zarr) : zarr(zarr)
    {
      if (zarr.get_name() != "time")
      {
        const std::string errmsg = "name of storage meant for "
                                   "timeobserver is not called 'time'";
        throw std::invalid_argument(errmsg);
      }
    }

    void observe_state(const size_t ngbxs,
                       const Kokkos::View<GridBox *> h_gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordinateStorage instance */
    {
      const auto &gbx = h_gridboxes(0);
      zarr.value_to_storage(gbx.state.time);
    }
};

class GridBoxIndexObserver
{
  private:
    CoordinateStorage<unsigned int> &zarr;
  
  public:
    GridBoxIndexObserver(CoordinateStorage<unsigned int> &zarr) : zarr(zarr)
    {
      if (zarr.get_name() != "gbxindex")
      {
        const std::string errmsg = "name of storage meant for gridbox index"
                                   " observer is not called 'gbxindex'";
        throw std::invalid_argument(errmsg);
      }
    }

    void observe_state(const size_t ngbxs,
                       const Kokkos::View<GridBox *> h_gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordinateStorage instance */
    {
      if (zarr.get_ndata() == 0)
      {
        for (size_t ii(0); ii<ngbxs; ++ii)
        { 
          zarr.value_to_storage(h_gridboxes(ii).gbxindex);
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

    void observe_state(const size_t ngbxs,
                       const Kokkos::View<GridBox *> h_gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordinateStorage instance */
    {
      for (size_t ii(0); ii<ngbxs; ++ii)
      {
        size_t nsupers = h_gridboxes(ii).span4SDsinGBx.size();
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

    void observe_state(const size_t ngbxs,
                       const Kokkos::View<GridBox *> h_gridboxes) const
    /* observe time of 0th gridbox and write it to an array
    as determined by the CoordinateStorage instance */
    {
      for (size_t ii(0); ii<ngbxs; ++ii)
      {
        const double moment = massnthmoment(h_gridboxes(ii).span4SDsinGBx,
                                            nth_moment);
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

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  /* observe superdroplets by writing their data to contigious
  ragged represented arrays as determined by the 
  ContiguousRaggedSDStorage instance */
  {
    size_t nsupers(0);
    for (size_t ii(0); ii<ngbxs; ++ii)
    {
      for (auto &SDinGBx : h_gridboxes(ii).span4SDsinGBx)
      {
        zarr.data_to_contigraggedarray(SDinGBx.superdrop);
        ++nsupers;
      }
    }
    zarr.contigraggedarray_count(nsupers);
  }
};

class SDsGbxindexObserver
{
private:
  ContiguousRaggedSDStorage<SdgbxIntoStore> &zarr;

public:
  SDsGbxindexObserver(ContiguousRaggedSDStorage<SdgbxIntoStore> &zarr) : zarr(zarr) {}

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  /* observe superdroplets by writing their data to contigious
  ragged represented arrays as determined by the 
  ContiguousRaggedSDStorage instance */
  {
    size_t nsupers(0);
    for (size_t ii(0); ii<ngbxs; ++ii)
    {
      for (auto &SDinGBx : h_gridboxes(ii).span4SDsinGBx)
      {
        zarr.data_to_contigraggedarray<unsigned int>(SDinGBx.sd_gbxindex);
        ++nsupers;
      }
    }
    zarr.contigraggedarray_count(nsupers);
  }
};

#endif // INTOSTORE_OBSERVERS_HPP