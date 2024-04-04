/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: massmoments_observer.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 4th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output the mass moments of the droplet size distribution in each gridbox
 * to individual arrays in a dataset a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./collect_data_for_dataset.hpp"
#include "./observers.hpp"
#include "./write_to_dataset_observer.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr2/dataset.hpp"

/* struct satifying CollectDataForDataset for collecting the 0th, 1st and 2nd mass moments of the
 * (rain)droplet distirubtion in each gridbox */
template <typename Store>
struct CollectMassMoments {
  // TODO(CB)

 public:
  struct Functor {
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {}  // TODO(CB)
  };

  Functor get_functor(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    return Functor{};  // TODO(CB)
  }

  void write_to_arrays(const Dataset<Store> &dataset) const {}  // TODO(CB)

  void write_to_ragged_arrays(const Dataset<Store> &dataset) const {}

  void write_arrayshapes(const Dataset<Store> &dataset) const {}  // TODO(CB)

  void write_ragged_arrayshapes(const Dataset<Store> &dataset) const {}

  void reallocate_views(const size_t sz) const {}
};

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
template <typename Store>
inline Observer auto MassMomentsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                         const int maxchunk, const size_t ngbxs) {
  const auto massmoments = 1;
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesTeamPolicyFunc{}, dataset, massmoments);
  return WriteToDatasetObserver(interval, parallel_write);
}

/* constructs observer which writes mass moments of raindroplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the WritetoDatasetObserver class */
template <typename Store>
inline Observer auto MassMomentsRaindropsObserver(const unsigned int interval,
                                                  const Dataset<Store> &dataset, const int maxchunk,
                                                  const size_t ngbxs) {
  const auto massmoments_raindrops = 1;
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesTeamPolicyFunc{}, dataset, massmoments_raindrops);
  return WriteToDatasetObserver(interval, parallel_write);
}

#endif  // LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_
