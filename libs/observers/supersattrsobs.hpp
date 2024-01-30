/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: supersattrsobs.hpp
 * Project: observers
 * Created Date: Monday 23rd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 24th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 */

#ifndef LIBS_OBSERVERS_SUPERSATTRSOBS_HPP_
#define LIBS_OBSERVERS_SUPERSATTRSOBS_HPP_

#include <concepts>
#include <iostream>
#include <memory>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr/contigraggedstorage.hpp"
#include "zarr/superdropattrsbuffers.hpp"
#include "zarr/superdropsbuffers.hpp"

/* constructs observer of the attributes of
all superdroplets in each gridbox with a
constant timestep 'interval' using an instance
of the DoStateObs class */
template <SuperdropsBuffers Buffers>
inline Observer auto SupersAttrsObserver(const unsigned int interval, FSStore &store,
                                         const int maxchunk, Buffers buffers);

/* observe superdroplets by writing their (attributes')
data to contigious ragged represented arrays as
determined by the ContigRaggedStorage instance */
template <SuperdropsBuffers Buffers>
struct DoSupersAttrsObs {
 private:
  using store_type = ContigRaggedStorage<Buffers>;
  std::shared_ptr<store_type> zarr;

 public:
  DoSupersAttrsObs(FSStore &store, const int maxchunk, Buffers buffers)
      : zarr(std::make_shared<store_type>(store, maxchunk, buffers)) {}

  void before_timestepping(const viewh_constgbx h_gbxs) const {
    std::cout << "observer includes SupersAttrsObserver\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewh_constgbx h_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(totsupers);
  }

  /* writes some variables from gridbox state
  to 2-D zarr storages as determined by the
  StateBuffers struct */
  void at_start_step(const viewd_constsupers d_totsupers) const {
    auto h_totsupers = Kokkos::create_mirror_view(
        d_totsupers);  // mirror of supers in case view is on device memory
    Kokkos::deep_copy(h_totsupers, d_totsupers);

    const size_t totnsupers(h_totsupers.extent(0));

    size_t obs_nsupers(0);
    for (size_t kk(0); kk < totnsupers; ++kk) {
      zarr->data_to_raggedstorage(h_totsupers(kk));
      ++obs_nsupers;
    }

    zarr->raggedarray_count(obs_nsupers);
  }

  void at_start_step(const unsigned int t_mdl, const Gridbox &gbx) const {}
};

/* constructs observer of the attributes of
all superdroplets in each gridbox with a
constant timestep 'interval' using an instance
of the DoStateObs class */
template <SuperdropsBuffers Buffers>
inline Observer auto SupersAttrsObserver(const unsigned int interval, FSStore &store,
                                         const int maxchunk, Buffers buffers) {
  const auto obs = DoSupersAttrsObs(store, maxchunk, buffers);
  return ConstTstepObserver(interval, obs);
}

#endif  // LIBS_OBSERVERS_SUPERSATTRSOBS_HPP_
