/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
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
 * File Description:
 * Observers to output nsupers (per gridbox
 * or total in domain) to array in a zarr
 * file system storage
 */

#ifndef LIBS_OBSERVERS_NSUPERSOBS_HPP_
#define LIBS_OBSERVERS_NSUPERSOBS_HPP_

#include <concepts>
#include <iostream>
#include <memory>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/supersingbx.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr/onedstorage.hpp"
#include "zarr/twodstorage.hpp"

namespace dlc = dimless_constants;

/* constructs observer of nsupers in each gridbox
with a constant timestep 'interval' using an
instance of the DoNsupersObs class */
inline Observer auto NsupersObserver(const unsigned int interval, FSStore &store,
                                     const int maxchunk, const size_t ngbxs);

/* constructs observer of nsupers that are raindrops
(r > rlim) in each gridbox, with a constant timestep
'interval' using an instance of the DoNrainsupersObs
class */
inline Observer auto NrainsupersObserver(const unsigned int interval, FSStore &store,
                                         const int maxchunk, const size_t ngbxs);

/* constructs observer of nsupers in entire domain,
with a constant timestep 'interval' using an
instance of the DoTotNsupersObs class */
inline Observer auto TotNsupersObserver(const unsigned int interval, FSStore &store,
                                        const int maxchunk, const size_t ngbxs);

/* observe nsupers in each gridbox and write
it to an array 'zarr' store as determined by
the 2DStorage instance */
class DoNsupersObs {
 private:
  using store_type = TwoDStorage<size_t>;
  std::shared_ptr<store_type> zarr;

 public:
  DoNsupersObs(FSStore &store, const int maxchunk, const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk, "nsupers", "<u8", " ", 1, "gbxindex",
                                          ngbxs)) {
    zarr->is_name("nsupers");
    zarr->is_dim1(ngbxs, "gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const {
    std::cout << "observer includes NsupersObserver\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {}

  /* gets number of superdrops for a gridbox
  and writes it to a 2-D zarr storage */
  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {
    const size_t nsupers(gbx.supersingbx.nsupers());
    zarr->value_to_storage(nsupers);
  }
};

/* constructs observer of nsupers in each gridbox
with a constant timestep 'interval' using an
instance of the DoNSupersObs class */
inline Observer auto NsupersObserver(const unsigned int interval, FSStore &store,
                                     const int maxchunk, const size_t ngbxs) {
  const auto obs = DoNsupersObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

/* observation for nsupers that are
raindrops (r > rlim) in each gridbox
which writes it to an array 'zarr' store
as determined by the 2DStorage instance */
class DoNrainsupersObs {
 private:
  using store_type = TwoDStorage<size_t>;
  std::shared_ptr<store_type> zarr;

  void nrainsupers_to_storage(const Gridbox &gbx) const;

 public:
  DoNrainsupersObs(FSStore &store, const int maxchunk, const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk, "nrainsupers", "<u8", " ", 1, "gbxindex",
                                          ngbxs)) {
    zarr->is_name("nrainsupers");
    zarr->is_dim1(ngbxs, "gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const {
    std::cout << "observer includes NrainsupersObserver\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {}

  /* Counts number of "raindrop-like" superdrops in
  gridbox and writes total number to 2-D zarr storage */
  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {
    nrainsupers_to_storage(gbx);
  }
};

/* constructs observer of nsupers that are raindrops
(r > rlim) in each gridbox, with a constant timestep
'interval' using an instance of the DoNrainsupersObs
class */
inline Observer auto NrainsupersObserver(const unsigned int interval, FSStore &store,
                                         const int maxchunk, const size_t ngbxs) {
  const auto obs = DoNrainsupersObs(store, maxchunk, ngbxs);
  return ConstTstepObserver(interval, obs);
}

/* observation of total nsupers in domain
(same as raggedcount) which writes it to an
array 'zarr' store as determined by
the 1DStorage instance */
class DoTotNsupersObs {
 private:
  using store_type = OneDStorage<size_t>;
  std::shared_ptr<store_type> zarr;

 public:
  DoTotNsupersObs(FSStore &store, const int maxchunk)
      : zarr(std::make_shared<store_type>(store, maxchunk, "totnsupers", "<u8", " ", 1)) {
    zarr->is_name("totnsupers");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const {
    std::cout << "observer includes TotNsupersObserver\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(totsupers);
  }

  /* gets number of superdrops in domain (from 0th gridbox
  metadata on superdrops) and writes it to a 1-D zarr storage */
  void at_start_step(const viewd_constsupers totsupers) const {
    const size_t totnsupers(totsupers.extent(0));
    zarr->value_to_storage(totnsupers);
  }

  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {}
};

/* constructs observer of total nsupers in domain
with a constant timestep 'interval' using an
instance of the DoTotNsupersObs class */
inline Observer auto TotNsupersObserver(const unsigned int interval, FSStore &store,
                                        const int maxchunk) {
  const auto obs = DoTotNsupersObs(store, maxchunk);
  return ConstTstepObserver(interval, obs);
}

#endif  // LIBS_OBSERVERS_NSUPERSOBS_HPP_
