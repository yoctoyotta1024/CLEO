/*
 * ----- CLEO -----
 * File: nsupersobs.hpp
 * Project: observers
 * Created Date: Friday 20th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 22nd October 2023
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
#include "gridboxes/gridbox.hpp"
#include "zarr/twodstorage.hpp"
#include "zarr/onedstorage.hpp"

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
the 2DVarStorage instance */
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
                     const viewh_constgbx h_gbxs) const
  /* converts integer model timestep to dimensionless time,
  then writes to zarr coordinate storage */
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
as determined by the 2DVarStorage instance */
{
private:
  using store_type = TwoDStorage<size_t>;
  std::shared_ptr<store_type> zarr;
  
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
                     const viewh_constgbx h_gbxs) const
  /* converts integer model timestep to dimensionless time,
  then writes to zarr coordinate storage */
  {
    constexpr double rlim(40e-6 / dlc::R0); // dimless minimum radius of raindrop

    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      auto d_supers = h_gbxs(ii).supersingbx();
      auto h_supers = Kokkos::create_mirror_view(d_supers); // mirror of supers in gridbox in case view is on device memory
      Kokkos::deep_copy(h_supers, d_supers);

      size_t nrainsupers(0);
      for (size_t kk(0); kk < h_supers.extent(0); ++kk)
      {
        const auto superdrop = h_supers(kk);
        if (superdrop.get_radius() >= rlim)
        {
          ++nrainsupers;
        }
      }
      zarr->value_to_storage(nrainsupers);
    }
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
the 2DVarStorage instance */
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
                     const viewh_constgbx h_gbxs) const
  /* converts integer model timestep to dimensionless time,
  then writes to zarr coordinate storage */
  {
    const auto gbx0 = h_gbxs(0);
    const size_t totnsupers(gbx0.supersingbx.domaintotnsupers());

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