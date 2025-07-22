/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: gbxindex_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output gridbox indexes at the start of each simulation to an array in a dataset
 */

#ifndef LIBS_OBSERVERS_GBXINDEX_OBSERVER_HPP_
#define LIBS_OBSERVERS_GBXINDEX_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "observers/observers.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "zarr/buffer.hpp"
#include "zarr/xarray_zarr_array.hpp"

/**
 * @struct GbxIndexFunctor
 * @brief Functor for copying gridbox indexes to a view in device memory.
 */
struct GbxIndexFunctor {
  viewd_constgbx d_gbxs;                       /**< View of gridboxes. */
  Buffer<uint32_t>::mirrorviewd_buffer d_data; /**< Mirror view on device for gridbox indexes. */

  /**
   * @brief Constructor for GbxIndexFunctor.
   * @param d_gbxs View of gridboxes on device.
   * @param d_data Mirror view on device of memory to copy gridbox indexes into.
   */
  GbxIndexFunctor(const viewd_constgbx d_gbxs, Buffer<uint32_t>::mirrorviewd_buffer d_data)
      : d_gbxs(d_gbxs), d_data(d_data) {}

  /**
   * @brief Operator is functor for within a Kokkos::parallel_for loop.
   *
   * Copies a gridbox's gbxindex to d_data view (for each gridbox in parallel).
   *
   * @param ii index for the gridbox in the d_gbxs and d_data views.
   */
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii) const { d_data(ii) = d_gbxs(ii).get_gbxindex(); }
};

/* @class GbxindexObserver
 * @brief Observer to output gridbox indexes to a 1-D array as coordinate of an xarray dataset.
 * @tparam Dataset Type of dataset.
 * @tparam Store Type of store for the dataset.
 */
template <typename Dataset, typename Store>
class GbxindexObserver {
 private:
  Dataset &dataset; /**< Dataset to write gridbox index data to. */
  std::shared_ptr<XarrayZarrArray<Store, uint32_t>>
      xzarr_ptr; /**< Pointer to gridbox index array in dataset. */

  /**
   * @brief Collects gridbox indexes from g_gbxs into a host memory view.
   * @param d_gbxs View of the gridboxes in device memory.
   * @return View in host memory of gridbox index of every gridbox in d_gbxs.
   */
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
  /**
   * @brief Constructor for GbxindexObserver.
   * @param dataset Dataset to write gridbox index data to.
   * @param store Store which dataset writes to.
   * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
   * @param ngbxs Number of gridboxes in final array.
   */
  GbxindexObserver(Dataset &dataset, Store &store, const size_t maxchunk, const size_t ngbxs)
      : dataset(dataset),
        xzarr_ptr(std::make_shared<XarrayZarrArray<Store, uint32_t>>(
            dataset.template create_coordinate_array<uint32_t>("gbxindex", "", 1, maxchunk,
                                                               ngbxs))) {}

  ~GbxindexObserver() { dataset.write_arrayshape(xzarr_ptr); }

  /**
   * @brief Obsevers the gridboxes' gbxindexes before timestepping.
   *
   * Write the gbxindex of every gridbox in d_gbxs to the gbxindex array in the dataset and assert
   * the size of the gbxindex dimension in the dataset is correct.
   *
   * @param d_gbxs View of gridboxes on device.
   */
  void before_timestepping(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    std::cout << "observer includes gbxindex observer\n";

    auto h_data = collect_gbxindexes(d_gbxs);
    dataset.write_to_array(xzarr_ptr, h_data);
    // assert((dataset.get_dimension("gbxindex") == h_data.extent(0)) &&
    //        "inconsistent size of gbxindex data and dataset dimension");
  }

  /**
   * @brief Placeholder for after timestepping functionality and to make class satisfy observer
   * concept.
   */
  void after_timestepping() const {}

  /**
   * @brief Placeholder for functionality at the start of each timestep and to make class satisfy
   * observer concept.
   * @param t_mdl Current model timestep.
   * @param d_gbxs View of gridboxes on device.
   * @param d_supers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const subviewd_constsupers d_supers) const {}

  /**
   * @brief Get null monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer that does nothing
   */
  SDMMonitor auto get_sdmmonitor() const { return NullSDMMonitor{}; }

  /**
   * @brief Returns the timestep of the next observation.
   *
   * No observation during timestepping so function returns the largest possible timestep (largest
   * unsigned integer).
   *
   * @param t_mdl Current model timestep.
   * @return Next observation timestep.
   */
  unsigned int next_obs(const unsigned int t_mdl) const { return LIMITVALUES::uintmax; }

  /**
   * @brief Checks if the current timestep is an observation timestep.
   *
   * No observation during timestepping so function always returns false.
   *
   * @param t_mdl Current model timestep.
   * @return boolean, false.
   */
  bool on_step(const unsigned int t_mdl) const { return false; }
};

#endif  // LIBS_OBSERVERS_GBXINDEX_OBSERVER_HPP_
