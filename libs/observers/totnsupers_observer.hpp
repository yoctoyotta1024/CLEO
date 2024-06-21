/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: totnsupers_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 21st June 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output total number of superdroplets at the start of each timestep to an array in a
 * dataset
 */

#ifndef LIBS_OBSERVERS_TOTNSUPERS_OBSERVER_HPP_
#define LIBS_OBSERVERS_TOTNSUPERS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "observers/consttstep_observer.hpp"
#include "observers/observers.hpp"
#include "superdrops/sdmmonitor.hpp"
#include "zarr/buffer.hpp"
#include "zarr/dataset.hpp"
#include "zarr/xarray_zarr_array.hpp"

/**
 * @class DoTotNsupersObs
 * @brief Template class for functionality to observe the total number of superdroplets at the start
 * of each timestep and write it to a Zarr array in an Xarray dataset.
 * @tparam Store Type of store for dataset.
 */
template <typename Store>
class DoTotNsupersObs {
 private:
  Dataset<Store> &dataset; /**< dataset to write totnsupers data to */
  std::shared_ptr<XarrayZarrArray<Store, uint32_t>> xzarr_ptr; /**< pointer to totnsupers array */

  /**
   * @brief Write out the total number of superdroplets in totsupers view at the start of a timestep
   * to an array in the dataset.
   *
   * _Note:_ conversion of totnsupers from size_t (arch dependent usually 8 bytes) to shorter 4
   * byte, unsigned int (unit32_t).
   *
   * @param totsupers View of all the superdroplets to count (on device memory but metadata for
   * extent of view is available on host).
   */
  void at_start_step(const viewd_constsupers totsupers) const {
    const auto data = static_cast<uint32_t>(totsupers.extent(0));
    dataset.write_to_array(xzarr_ptr, data);
  }

 public:
  /**
   * @brief Constructor for DoTotNsupersObs.
   * @param dataset Dataset to write totnsupers data to.
   * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
   */
  DoTotNsupersObs(Dataset<Store> &dataset, const size_t maxchunk)
      : dataset(dataset),
        xzarr_ptr(std::make_shared<XarrayZarrArray<Store, uint32_t>>(
            dataset.template create_array<uint32_t>("totnsupers", "", 1, {maxchunk}, {"time"}))) {}

  /**
   * @brief Destructor for DoTotNsupersObs.
   */
  ~DoTotNsupersObs() { dataset.write_arrayshape(xzarr_ptr); }

  /**
   * @brief Placeholder for before timestepping functionality and to make class satisfy observer
   * concept.
   */
  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes totnsupers observer\n";
  }

  /**
   * @brief Placeholder for after timestepping functionality and to make class satisfy observer
   * concept.
   */
  void after_timestepping() const {}

  /**
   * @brief Adapter to call at start step function which writes the total number of superdroplets in
   * totsupers view to the array in the dataset.
   *
   * @param t_mdl Current model timestep.
   * @param d_gbxs View of gridboxes on device.
   * @param totsupers View of superdrops on device.
   */
  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(totsupers);
  }

  /**
   * @brief Get null monitor for SDM processes from observer.
   *
   * @return monitor 'mo' of the observer that does nothing
   */
  SDMMonitor auto get_sdmmonitor() const { return NullSDMMonitor{}; }
};

/**
 * @brief Constructs an observer which writes the total number of superdroplets at start of each
 * observation timestep to a 1-D array with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @return Constructed type satisfying observer concept.
 */
template <typename Store>
inline Observer auto TotNsupersObserver(const unsigned int interval, Dataset<Store> &dataset,
                                        const size_t maxchunk) {
  return ConstTstepObserver(interval, DoTotNsupersObs(dataset, maxchunk));
}

#endif  // LIBS_OBSERVERS_TOTNSUPERS_OBSERVER_HPP_
