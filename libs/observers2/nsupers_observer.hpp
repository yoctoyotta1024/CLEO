/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: nsupers_observer.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 2nd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output  at the start of
 * each timestep to an array in a dataset
 */

#ifndef LIBS_OBSERVERS2_NSUPERS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_NSUPERS_OBSERVER_HPP_

#include <concepts>

#include "./do_write_gridboxes.hpp"
#include "./generic_write_gridbox_to_array.hpp"
#include "./observers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/buffer.hpp"
#include "zarr2/dataset.hpp"

/* Operator is functor to perform copy of number of superdrops in each gridbox to d_data
in parallel. Note conversion of nsupers from size_t (8 bytes) to single precision (4 bytes
unsigned integer) in output */
struct NsupersFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto nsupers = static_cast<uint32_t>(d_gbxs(ii).supersingbx.nsupers());
    d_data(ii) = nsupers;
  }
};

/* returns WriteGridboxToArray which writes the number of superdrops in each
gridbox to an array in a dataset in a store */
template <typename Store>
WriteGridboxToArray<Store> auto NsupersWriter(Dataset<Store> &dataset, const size_t maxchunk,
                                              const size_t ngbxs) {
  return GenericWriteGridboxToXarray<Store, uint32_t, NsupersFunc>(dataset, "nsupers", "", "<u4", 1,
                                                                   maxchunk, ngbxs, NsupersFunc{});
}

/* constructs observer which writes number of superdrops in each gridbox with a constant
timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto NsupersObserver(const unsigned int interval, Dataset<Store> &dataset,
                                     const int maxchunk, const size_t ngbxs) {
  const WriteGridboxToArray<Store> auto nsuperswriter = NsupersWriter(dataset, maxchunk, ngbxs);
  const auto looppolicy = ParallelGbxsRangePolicy{};
  return ConstTstepObserver(interval, DoWriteGridboxes(dataset, nsuperswriter, looppolicy));
}

#endif  // LIBS_OBSERVERS2_NSUPERS_OBSERVER_HPP_

// TODO(CB) totnsupers in domain observer
