// Author: Clara Bayley
// File: "gridboxes_intostore.hpp"
/* structures obeying the ObserveGBxs
concept for various ways of observing
a gridbox which ends up writing data
into a (zarr) store on disk */

#ifndef GRIDBOXES_INTOSTORE_HPP
#define GRIDBOXES_INTOSTORE_HPP

#include <string>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "./observegbxs.hpp"
#include "zarrstorage/thermostatestorage.hpp"
#include "zarrstorage/sdattributes_intostore.hpp"
#include "zarrstorage/contigraggedsdstorage.hpp"
#include "zarrstorage/singlevarstorage.hpp"
#include "zarrstorage/massmomentsstorage.hpp"
#include "sdmgridboxes/gridbox.hpp"

class ObserveThermoState
/* observe thermostate of each gridbox by
writing it to arrays in a zarr store as
determined by the ThermoStateStorage instance */
{
private:
  ThermoStateStorage &zarr;

public:
  ObserveThermoState(ThermoStateStorage &zarr) : zarr(zarr){}

  void prepare() const {}

  void operator()(const size_t ngbxs,
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
class ObserveSDsAttributes
/* observe superdroplets by writing their (attributes')
data to contigious ragged represented arrays as
determined by the ContiguousRaggedSDStorage instance */
{
private:
  ContiguousRaggedSDStorage &zarr;

public:
  ObserveSDsAttributes(ContiguousRaggedSDStorage &zarr) : zarr(zarr) {}

  void prepare() const {}

  void operator()(const size_t ngbxs,
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

class ObserveSDsGbxindex
/* observe gridbox index of each superdroplet and write to
zarr storage in a contigious ragged represented array as
determined by the ContiguousRaggedSDStorage instance */
{
private:
  ContiguousRaggedSDStorage<SdgbxIntoStore> &zarr;

public:
  ObserveSDsGbxindex(auto &zarr) : zarr(zarr) {}
  
  void prepare() const {}

  void operator()(const size_t ngbxs,
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

class ObserveTime
/* observe time of 0th gridbox and write it
to an array 'zarr' store as determined by
the CoordinateStorage instance */
{
private:
  CoordinateStorage<double> &zarr;

public:
  ObserveTime(CoordinateStorage<double> &zarr) : zarr(zarr)
  {
    zarr.is_name("time");
  }

  void prepare() const { zarr.is_name("time"); }

  void operator()(const size_t ngbxs,
                  const Kokkos::View<GridBox *> h_gridboxes) const

  {
    const auto &gbx = h_gridboxes(0);
    zarr.value_to_storage(gbx.state.time);
  }
};

class ObserveGridBoxIndex
/* observe the gbxindex of each gridbox and
write it to an array 'zarr' store as determined
by the CoordinateStorage instance */
{
private:
  CoordinateStorage<unsigned int> &zarr;

public:
  ObserveGridBoxIndex(CoordinateStorage<unsigned int> &zarr)
      : zarr(zarr)
  {
    zarr.is_name("gbxindex");
  }

  void prepare() const { zarr.is_name("gbxindex"); }

  void operator()(const size_t ngbxs,
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

class ObserveNsupersPerGridBox
/* observe number of superdroplets in each gridbox
and write to 'zarr', a 2D array in a zarr store */
{
private:
  TwoDStorage<size_t> &zarr;

public:
  ObserveNsupersPerGridBox(TwoDStorage<size_t> &zarr,
                           const size_t ngbxs)
      : zarr(zarr)
  {
    zarr.is_name("nsupers");
    zarr.is_dim1(ngbxs, "gbxindex");
  }

  void prepare() const { zarr.is_name("nsupers"); }

  void operator()(const size_t ngbxs,
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

class ObserveNthMassMoment
/* observe nth mass moment of (real) droplets
distribution in each gridbox through by
writing data from 'massmoment' function
to an array in a zarr storage 'zarr' */
{
private:
  const int nth_moment;
  TwoDStorage<double> &zarr;

public:
  ObserveNthMassMoment(TwoDStorage<double> &zarr,
                        const int nth_moment,
                        const size_t ngbxs)
      : nth_moment(nth_moment),
        zarr(zarr)
  {
    const std::string name("mom" + std::to_string(nth_moment));
    zarr.is_name(name);
    zarr.is_dim1(ngbxs, "gbxindex");
  }

  void prepare() const
  {
    const std::string name("mom" + std::to_string(nth_moment));
    zarr.is_name(name);
  }

  void operator()(const size_t ngbxs,
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

class ObserveNthRainMassMoment
/* observe nth mass moment of raindroplets
distribution in each gridbox through by
writing data from 'rainmassmoment' function
to an array in a zarr storage 'zarr' */
{
private:
  const int nth_moment;
  TwoDStorage<double> &zarr;

public:
  ObserveNthRainMassMoment(TwoDStorage<double> &zarr,
                            const int nth_moment,
                            const size_t ngbxs)
      : nth_moment(nth_moment),
        zarr(zarr)
  {
    const std::string name("rainmom" + std::to_string(nth_moment));
    zarr.is_name(name);
    zarr.is_dim1(ngbxs, "gbxindex");
  }

  void prepare() const
  {
    const std::string name("rainmom" + std::to_string(nth_moment));
    zarr.is_name(name);
  } 

  void operator()(const size_t ngbxs,
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

class ObserveNRainsupersPerGridBox
/* observe number of superdroplets in each gridbox
and write to 'zarr', a 2D array in a zarr store */
{
private:
  TwoDStorage<size_t> &zarr;

public:
  ObserveNRainsupersPerGridBox(TwoDStorage<size_t> &zarr,
                           const size_t ngbxs)
      : zarr(zarr)
  {
    zarr.is_name("nrainsupers");
    zarr.is_dim1(ngbxs, "gbxindex");
  }

  void prepare() const { zarr.is_name("nsupers"); }

  void operator()(const size_t ngbxs,
                  const Kokkos::View<GridBox *> h_gridboxes) const
  {
    constexpr double rlim(40e-6/dlc::R0); // dimless minimum radius of precip
    
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      size_t nrainsupers(0);
      for (const auto &SDinGBx : h_gridboxes(ii).span4SDsinGBx)
      {
        if (SDinGBx.superdrop.radius >= rlim)
        {
          ++nrainsupers;
        }
      }
      zarr.value_to_storage(nrainsupers);
    }
    ++zarr.nobs;
  }
};

#endif // GRIDBOXES_INTOSTORE_HPP