/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: gbxindex_observer.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 28th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output gridbox indexes at the start of each simulation to an array in a dataset
 */

#ifndef LIBS_OBSERVERS2_GBXINDEX_OBSERVER_HPP_
#define LIBS_OBSERVERS2_GBXINDEX_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <functional>
#include <iostream>
#include <memory>
#include <utility>

#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"

struct GbxIndexFunctor {
  viewd_constgbx d_gbxs;  // view of gridboxes
  Buffer<uint32_t>::mirrorviewd_buffer
      d_data;  // mirror view on device for gbxindex of every gridbox

  GbxIndexFunctor(const viewd_constgbx d_gbxs, Buffer<uint32_t>::mirrorviewd_buffer d_data)
      : d_gbxs(d_gbxs), d_data(d_data) {}

  // Functor operator to perform copy of gbxindex of each gridbox to d_data in parallel
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii) const { d_data(ii) = d_gbxs(ii).get_gbxindex(); }
};

/* template observer which writes gbxindex for every gridbox out to a 1-D array
as a coordinate of an xarray dataset */
template <typename Store>
class GbxindexObserver {
 private:
  Dataset<Store> &dataset;  ///< dataset to write gbxindex data to
  std::shared_ptr<XarrayZarrArray<Store, uint32_t>> xzarr_ptr;  ///< pointer to gbxindex array

  /* returns a view in the host memor of the gbxindex of every gridbox in d_gbxs */
  Buffer<uint32_t>::viewh_buffer collect_gbxindexes(const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    auto h_data = Buffer<uint32_t>::viewh_buffer("h_data", ngbxs);
    auto d_data = Kokkos::create_mirror_view(ExecSpace(), h_data);
    Kokkos::parallel_for("collect_gbxs_data", Kokkos::RangePolicy<ExecSpace>(0, ngbxs),
                         GbxIndexFunctor(d_gbxs, d_data));
    Kokkos::deep_copy(h_data, d_data);
    return h_data;
  }

 public:
  GbxindexObserver(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs)
      : dataset(dataset),
        xzarr_ptr(std::make_shared<XarrayZarrArray<Store, uint32_t>>(
            dataset.template create_coordinate_array<uint32_t>("gbxindex", "", "<u4", 1, maxchunk,
                                                               ngbxs))) {}

  ~GbxindexObserver() { dataset.write_arrayshape(xzarr_ptr); }

  /* write the gbxindex of every gridbox in d_gbxs to the gbxindex array in the dataset and assert
  the size of the gbxindex dimension in the dataset is correct */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes gbxindex observer\n";
    auto h_data = collect_gbxindexes(d_gbxs);
    dataset.write_to_array(xzarr_ptr, h_data);
    assert((dataset.get_dimension("gbxindex") == h_data.extent(0)) &&
           "inconsistent size of gbxindex data and dataset dimension");
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {}

  unsigned int next_obs(const unsigned int t_mdl) const { return LIMITVALUES::uintmax; }

  bool on_step(const unsigned int t_mdl) const { return false; }
};

#endif  // LIBS_OBSERVERS2_GBXINDEX_OBSERVER_HPP_
