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
 * Last Modified: Saturday 30th March 2024
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

#include "./observers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "./write_gridboxes.hpp"
#include "zarr2/dataset.hpp"

/* constructs observer which writes number of superdrops in each gridbox with a constant
timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto NsupersObserver(const unsigned int interval, Dataset<Store> &dataset,
                                     const int maxchunk, const size_t ngbxs) {
  const WriteGridboxToArray<Store> auto nsuperswriter = NsupersWriter(dataset, maxchunk, ngbxs);
  return ConstTstepObserver(interval, WriteGridboxes(dataset, nsuperswriter));
}

#endif  // LIBS_OBSERVERS2_NSUPERS_OBSERVER_HPP_
