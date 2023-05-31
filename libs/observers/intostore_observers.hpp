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
    const std::string errmsg("name of storage is called " +
                             zarrname + ", but should be " + name);
    throw std::invalid_argument(errmsg);
  }
}

class ThermoStateObserver
/* observe thermostate of each gridbox by
writing it to arrays in a zarr store as
determined by the ThermoStateStorage instance */
{
private:
  ThermoStateStorage &zarr;

public:
  ConstIntervalStep on_step;

  ThermoStateObserver(const int obsstep, ThermoStateStorage &zarr)
      : zarr(zarr), on_step(obsstep) {}

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  {
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      zarr.thermodata_to_storage(h_gridboxes(ii).state);
    }
    ++zarr.nobs;
  }
};

template <typename ContiguousRaggedSDStorage>
class SDsAttributeObserver
/* observe superdroplets by writing their (attributes')
data to contigious ragged represented arrays as
determined by the ContiguousRaggedSDStorage instance */
{
private:
  ContiguousRaggedSDStorage &zarr;

public:
  ConstIntervalStep on_step;

  SDsAttributeObserver(const int obsstep,
                       ContiguousRaggedSDStorage &zarr)
      : zarr(zarr), on_step(obsstep) {}

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  {
    size_t totnsupers(0);
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      for (auto &SDinGBx : h_gridboxes(ii).span4SDsinGBx)
      {
        zarr.data_to_raggedstorage(SDinGBx.superdrop);
        ++totnsupers;
      }
    }
    zarr.raggedarray_count(totnsupers);
  }
};

class SDsGbxindexObserver
/* observe gridbox index of each superdroplet and write to
zarr storage in a contigious ragged represented array as
determined by the ContiguousRaggedSDStorage instance */
{
private:
  ContiguousRaggedSDStorage<SdgbxIntoStore> &zarr;

public:
  ConstIntervalStep on_step;

  SDsGbxindexObserver(const int obsstep, auto &zarr)
      : zarr(zarr), on_step(obsstep) {}

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  {
    size_t totnsupers(0);
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      for (auto &SDinGBx : h_gridboxes(ii).span4SDsinGBx)
      {
        zarr.data_to_raggedstorage<unsigned int>(SDinGBx.sd_gbxindex);
        ++totnsupers;
      }
    }
    zarr.raggedarray_count(totnsupers);
  }
};

class TimeObserver
/* observe time of 0th gridbox and write it
to an array 'zarr' store as determined by
the CoordinateStorage instance */
{
private:
  CoordinateStorage<double> &zarr;

public:
  ConstIntervalStep on_step;

  TimeObserver(const int obsstep, CoordinateStorage<double> &zarr)
      : zarr(zarr), on_step(obsstep)
  {
    check_zarrname(zarr.get_name(), "time");
  }

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const

  {
    const auto &gbx = h_gridboxes(0);
    zarr.value_to_storage(gbx.state.time);
  }
};

class GridBoxIndexObserver
/* observe the gbxindex of each gridbox and
write it to an array 'zarr' store as determined
by the CoordinateStorage instance */
{
private:
  CoordinateStorage<unsigned int> &zarr;

public:
  ConstIntervalStep on_step;

  GridBoxIndexObserver(const int obsstep,
                       CoordinateStorage<unsigned int> &zarr)
      : zarr(zarr), on_step(obsstep)
  {
    check_zarrname(zarr.get_name(), "gbxindex");
  }

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
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
/* observe number of superdroplets in each gridbox
and write to 'zarr', a 2D array in a zarr store */
{
private:
  TwoDStorage<size_t> &zarr;

public:
  ConstIntervalStep on_step;

  NsupersPerGridBoxObserver(const int obsstep,
                            TwoDStorage<size_t> &zarr)
      : zarr(zarr), on_step(obsstep)
  {
    check_zarrname(zarr.get_name(), "nsupers");
  }

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
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
/* observe nth mass moment of (real) droplets
distribution in each gridbox through by
writing data from 'massmoment' function
to an array in a zarr storage 'zarr' */
{
private:
  const int nth_moment;
  TwoDStorage<double> &zarr;

public:
  ConstIntervalStep on_step;

  NthMassMomentObserver(const int obsstep,
                        TwoDStorage<double> &zarr,
                        const int nth_moment)
      : nth_moment(nth_moment),
        zarr(zarr),
        on_step(obsstep)
  {
    const std::string name("mom" + std::to_string(nth_moment));
    check_zarrname(zarr.get_name(), name);
  }

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  {
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      zarr.value_to_storage(
          massmoment(h_gridboxes(ii).span4SDsinGBx, nth_moment));
    }
    ++zarr.nobs;
  }
};

class NthRainMassMomentObserver
/* observe nth mass moment of raindroplets
distribution in each gridbox through by
writing data from 'rainmassmoment' function
to an array in a zarr storage 'zarr' */
{
private:
  const int nth_moment;
  TwoDStorage<double> &zarr;

public:
  ConstIntervalStep on_step;

  NthRainMassMomentObserver(const int obsstep,
                            TwoDStorage<double> &zarr,
                            const int nth_moment)
      : nth_moment(nth_moment),
        zarr(zarr),
        on_step(obsstep)
  {
    const std::string name("rainmom" + std::to_string(nth_moment));
    check_zarrname(zarr.get_name(), name);
  }

  int get_interval() const { return on_step.get_interval(); }

  void observe_gridboxes(const size_t ngbxs,
                     const Kokkos::View<GridBox *> h_gridboxes) const
  {
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      zarr.value_to_storage(
          rainmassmoment(h_gridboxes(ii).span4SDsinGBx, nth_moment));
    }
    ++zarr.nobs;
  }
};

// class SurfacePrecipObserver
// /* observe surface precipitation and write it
// to an array in a zarr storage 'zarr' */
// {
// private:
//   TwoDStorage<double> &zarr;

//   const double zlim = 5.0 / dlc::COORD0; // dimless maximum z coord of precip
//   std::vector<size_t> surface_gbxindexes; //indexes of gbxs to include as surface

// public:
//   SurfacePrecipObserver(TwoDStorage<double> &zarr, gbxmaps)
//       : zarr(zarr)
//       {
//         surface_gbxindexes_pushback when gbxindex upper limit within zlim
//       }
//   {
//     check_zarrname(zarr.get_name(), "surfprecip");
//     TODO: need_to_set_units_and_sf_in_storage_zattrs_too
//   }

//   void observe_gridboxes(const size_t ngbxs,
//                      const Kokkos::View<GridBox *> h_gridboxes) const
//   {
//     for (size_t ii(0); ii < ngbxs; ++ii)
//     {
//       if (gbx.index not in list_of_gbxindexes_to_include zlim)
//       {
//         zarr.value_to_storage(0.0);
//       }
//       else
//       {
//         zarr.value_to_storage(surface_precipitation(h_gridboxes(ii), zlim));
//       }
//     }
//     ++zarr.nobs;
//   }
// };

#endif // INTOSTORE_OBSERVERS_HPP