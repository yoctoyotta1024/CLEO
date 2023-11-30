/*
 * ----- CLEO -----
 * File: nsupersobs.hpp
 * Project: observers
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 30th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Observers to output nsupers (per gridbox
 * or total in domain) to array in a zarr
 * file system storage
 */

#ifndef NSUPERSOBS_HPP 
#define NSUPERSOBS_HPP 

#include <concepts>
#include <memory>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "superdrops/superdrop.hpp"
#include "gridboxes/supersingbx.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/twodstorage.hpp"
#include "zarr/onedstorage.hpp"

namespace dlc = dimless_constants;

inline Observer auto
NsupersObserver(const unsigned int interval,
                FSStore &store,
                const int maxchunk, 
                const size_t ngbxs);
/* constructs observer of nsupers in each gridbox
with a constant timestep 'interval' using an
instance of the DoNsupersObs class */

inline Observer auto
NrainsupersObserver(const unsigned int interval,
                    FSStore &store,
                    const int maxchunk,
                    const size_t ngbxs);
/* constructs observer of nsupers that are raindrops
(r > rlim) in each gridbox, with a constant timestep
'interval' using an instance of the DoNrainsupersObs
class */

inline Observer auto
TotNsupersObserver(const unsigned int interval,
                   FSStore &store,
                   const int maxchunk,
                   const size_t ngbxs);
/* constructs observer of nsupers in entire domain,
with a constant timestep 'interval' using an
instance of the DoTotNsupersObs class */

class DoNsupersObs
/* observe nsupers in each gridbox and write
it to an array 'zarr' store as determined by
the 2DStorage instance */
{
private:
  using store_type = TwoDStorage<size_t>;
  std::shared_ptr<store_type> zarr;
  
public:
  DoNsupersObs(FSStore &store,
               const int maxchunk,
               const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "nsupers", "<u8",
                                          " ", 1,
                                          "gbxindex", ngbxs))
  {
    zarr->is_name("nsupers");
    zarr->is_dim1(ngbxs, "gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes NsupersObserver\n";
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const
  {
    at_start_step(h_gbxs);
  }

  void at_start_step(const viewh_constgbx h_gbxs) const
  /* gets number of superdrops for each gridbox and
  writes it to a 2-D zarr storage */
  {
    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      const size_t nsupers(h_gbxs(ii).supersingbx.nsupers());
      zarr->value_to_storage(nsupers);
    }
    ++(zarr->nobs);
  }
};

inline Observer auto
NsupersObserver(const unsigned int interval,
                FSStore &store,
                const int maxchunk,
                const size_t ngbxs)
/* constructs observer of nsupers in each gridbox
with a constant timestep 'interval' using an
instance of the DoNSupersObs class */
{
  const auto obs = DoNsupersObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

class DoNrainsupersObs
/* observation for nsupers that are
raindrops (r > rlim) in each gridbox
which writes it to an array 'zarr' store
as determined by the 2DStorage instance */
{
private:
  using store_type = TwoDStorage<size_t>;
  std::shared_ptr<store_type> zarr;
  
  void nrainsupers_to_storage(const viewh_constgbx h_gbxs) const;

public:
  DoNrainsupersObs(FSStore &store,
               const int maxchunk,
               const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "nrainsupers", "<u8",
                                          " ", 1,
                                          "gbxindex", ngbxs))
  {
    zarr->is_name("nrainsupers");
    zarr->is_dim1(ngbxs, "gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes NrainsupersObserver\n";
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const
  {
    at_start_step(h_gbxs);
  }

  void at_start_step(const viewh_constgbx h_gbxs) const
  /* Counts number of "raindrop-like" superdrops for each
  gridbox and writes total number to 2-D zarr storage */
  {
    nrainsupers_to_storage(h_gbxs);
    ++(zarr->nobs);
  }
};

inline Observer auto
NrainsupersObserver(const unsigned int interval,
                    FSStore &store,
                    const int maxchunk,
                    const size_t ngbxs)
/* constructs observer of nsupers that are raindrops
(r > rlim) in each gridbox, with a constant timestep
'interval' using an instance of the DoNrainsupersObs
class */
{
  const auto obs = DoNrainsupersObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

class DoTotNsupersObs
/* observation of total nsupers in domain
(same as raggedcount) which writes it to an
array 'zarr' store as determined by
the 1DStorage instance */
{
private:
  using store_type = OneDStorage<size_t>;
  std::shared_ptr<store_type> zarr;

public:
  DoTotNsupersObs(FSStore &store,
                  const int maxchunk)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                          "totnsupers", "<u8",
                                          " ", 1))
  {
    zarr->is_name("totnsupers");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes TotNsupersObserver\n";
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const
  {
    at_start_step(totsupers);
  }

  void at_start_step(const viewd_constsupers totsupers) const
  /* gets number of superdrops in domain (from 0th gridbox 
  metadata on superdrops) and writes it to a 1-D zarr storage */
  {
    const size_t totnsupers(totsupers.extent(0));

    zarr->value_to_storage(totnsupers);
  }
};


inline Observer auto
TotNsupersObserver(const unsigned int interval,
                    FSStore &store,
                    const int maxchunk)
/* constructs observer of total nsupers in domain
with a constant timestep 'interval' using an
instance of the DoTotNsupersObs class */
{
  const auto obs = DoTotNsupersObs(store, maxchunk);
  return ConstTstepObserver(interval, obs);
}

#endif // NSUPERSOBS_HPP 