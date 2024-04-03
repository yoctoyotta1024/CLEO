/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrops_observer.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 3rd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables related to Gridboxes' state at the start of
 * each timestep to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "./do_write_supers.hpp"
#include "./generic_write_supers_to_array.hpp"
#include "./observers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "zarr2/dataset.hpp"

/* Operator is functor to perform copy of xi for each superdroplet in totsupers view to d_data
in parallel. Note conversion of xi from size_t (arch dependent usually 8 bytes) to long
precision unsigned int (unit64_t) */
struct XiFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, const viewd_constsupers totsupers,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto xi = static_cast<uint32_t>(supers(kk).get_xi());
    d_data(kk) = xi;
  }
};

/* constructs observer which writes out variables from every superdroplet in the domain
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto SuperdropsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                        const int maxchunk) {
  const WriteGridboxToArray<Store> auto xi =
      GenericWriteSupersToXarray<Store, uint32_t, XiFunc, XarrayForSupersData>(
          dataset, "xi", "", "<u8", 1, maxchunk, XiFunc{});

  const auto obsfunc = DoWriteSupers(dataset, maxchunk, xi);
  return ConstTstepObserver(interval, obsfunc);
}

#endif  // LIBS_OBSERVERS2_SUPERDROPS_OBSERVER_HPP_
