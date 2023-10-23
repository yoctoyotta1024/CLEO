/*
 * ----- CLEO -----
 * File: massmomentsobs.hpp
 * Project: observers
 * Created Date: Sunday 22nd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 23rd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Observer to output nsupers per gridbox
 * to array in a zarr file system storage
 */

#ifndef MASSMOMENTSOBS_HPP
#define MASSMOMENTSOBS_HPP

#include <concepts>
#include <memory>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "superdrops/superdrop.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/massmomentsstorage.hpp"

inline Observer auto
MassMomentsObserver(const unsigned int interval,
                   FSStore &store,
                   const int maxchunk,
                   const size_t ngbxs);
/* constructs observer of the nth mass moment
'mom' in each gridbox with a constant
timestep 'interval' using an instance of the
DoMassMomentsObs class */

inline Observer auto
RainMassMomentsObserver(const unsigned int interval,
                        FSStore &store,
                        const int maxchunk,
                        const size_t ngbxs);
/* constructs observer of the nth mass moment
'mom' for raindrops (r>rlim) in each gridbox with a
constant timestep 'interval' using an instance
of the DoRainMassMomentsObs class */

class DoMassMomentsObs
/* observe nth mass moment in each gridbox and
write it to an array 'zarr' store as determined
by the 2DStorage instance */
{
private:
  using store_type = MassMomentsStorage<double>;
  std::shared_ptr<store_type> zarr;

  void massmoments(const mirrorh_constsupers h_supers) const;
  /* calculated 0th, 1st and 2nd moment of the (real) droplet mass
  distribution and then writes them to zarr storage. (I.e.
  0th, 3rd and 6th moment of the droplet radius distribution) */

public:
  DoMassMomentsObs(FSStore &store,
                   const int maxchunk,
                   const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "<f8", ngbxs, ""))
  {
    zarr->is_dim1(ngbxs, "gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes MassMomentsObserver\n";
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  /* deep copy if necessary (if superdrops are on device not
  host memory), then writes mass moments to 2-D zarr storages */
  {
    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      auto h_supers = h_gbxs(ii).supersingbx.hostcopy(); // deep copy if supers not on host memory
      massmoments(h_supers);
    }
    ++(zarr->nobs);
  }
};

inline Observer auto
MassMomentsObserver(const unsigned int interval,
                   FSStore &store,
                   const int maxchunk,
                   const size_t ngbxs)
/* constructs observer of the nth mass moment
'mom' in each gridbox with a constant
timestep 'interval' using an instance of the
DoMassMomentsObs class */
{
  const auto obs = DoMassMomentsObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

class DoRainMassMomentsObs
/* observe nth mass moment in each gridbox and
write it to an array 'zarr' store as determined
by the 2DStorage instance */
{
private:
  using store_type = MassMomentsStorage<double>;
  std::shared_ptr<store_type> zarr;

  void rainmassmoments(const mirrorh_constsupers h_supers) const;
  /* calculated 0th, 1st and 2nd moment of the (real) droplet mass
  distribution and then writes them to zarr storage. (I.e.
  0th, 3rd and 6th moment of the droplet radius distribution) */

public:
  DoRainMassMomentsObs(FSStore &store,
                   const int maxchunk,
                   const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "<f8", ngbxs,
                                          "rain"))
  {
    zarr->is_dim1(ngbxs, "gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes RainMassMomentsObserver\n";
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  /* deep copy if necessary (if superdrops are on device not
  host memory), then writes mass moments to 2-D zarr storages */
  {
    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      auto h_supers = h_gbxs(ii).supersingbx.hostcopy(); // deep copy if supers not on host memory
      rainmassmoments(h_supers);
    }
    ++(zarr->nobs);
  }
};

inline Observer auto
RainMassMomentsObserver(const unsigned int interval,
                        FSStore &store,
                        const int maxchunk,
                        const size_t ngbxs)
/* constructs observer of the nth mass moment
'mom' for raindrops (r>rlim) in each gridbox with a
constant timestep 'interval' using an instance
of the DoRainMassMomentsObs class */
{
  const auto obs = DoRainMassMomentsObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

#endif // MASSMOMENTSOBS_HPP
