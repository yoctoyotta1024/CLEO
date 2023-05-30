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
#include "./massmomentsstorage.hpp"
#include "sdmgridboxes/gridbox.hpp"

void check_zarrname(const std::string zarrname,
                    const std::string name)
{
  if (zarrname != name)
  {
    const std::string errmsg = "name of storage is called " +
                               zarrname + ", but should be " + name;
    throw std::invalid_argument(errmsg);
  }
}

class ThermoStateObserver
{
private:
  ThermoStateStorage &zarr;

public:
  ThermoStateObserver(ThermoStateStorage &zarr) : zarr(zarr) {}

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  /* observe thermostate by writing it to arrays
  as determined by the ThermoStateStorage instance */
  {
    for (size_t ii(0); ii < ngbxs; ++ii)
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
    check_zarrname(zarr.get_name(), "time");
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
    check_zarrname(zarr.get_name(), "gbxindex");
  }

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  /* observe time of 0th gridbox and write it to an array
  as determined by the CoordinateStorage instance */
  {
    if (zarr.get_ndata() == 0)
    {
      for (size_t ii(0); ii < ngbxs; ++ii)
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
    check_zarrname(zarr.get_name(), "nsupers");
  }

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  /* observe time of 0th gridbox and write it to an array
  as determined by the CoordinateStorage instance */
  {
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      size_t nsupers = h_gridboxes(ii).span4SDsinGBx.size();
      zarr.value_to_storage(nsupers);
    }

    ++zarr.nobs;
  }
};

class NthMassMomentObserver
{
private:
  const int nth_moment;
  TwoDStorage<double> &zarr;

public:
  NthMassMomentObserver(TwoDStorage<double> &zarr,
                          const int nth_moment)
      : nth_moment(nth_moment),
        zarr(zarr)
  {
    const std::string name("massmom" + std::to_string(nth_moment));
    check_zarrname(zarr.get_name(), name);
  }

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  /* observe time of 0th gridbox and write it to an array
  as determined by the CoordinateStorage instance */
  {
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      const double moment = massmoment(h_gridboxes(ii).span4SDsinGBx,
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
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      for (auto &SDinGBx : h_gridboxes(ii).span4SDsinGBx)
      {
        zarr.data_to_raggedstorage(SDinGBx.superdrop);
        ++nsupers;
      }
    }
    zarr.raggedarray_count(nsupers);
  }
};

class SDsGbxindexObserver
{
private:
  ContiguousRaggedSDStorage<SdgbxIntoStore> &zarr;

public:
  SDsGbxindexObserver(auto &zarr) : zarr(zarr) {}

  void observe_state(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  /* observe superdroplets by writing their data to contigious
  ragged represented arrays as determined by the
  ContiguousRaggedSDStorage instance */
  {
    size_t nsupers(0);
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      for (auto &SDinGBx : h_gridboxes(ii).span4SDsinGBx)
      {
        zarr.data_to_raggedstorage<unsigned int>(SDinGBx.sd_gbxindex);
        ++nsupers;
      }
    }
    zarr.raggedarray_count(nsupers);
  }
};

#endif // INTOSTORE_OBSERVERS_HPP