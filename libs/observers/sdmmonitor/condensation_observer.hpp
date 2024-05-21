/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: condensation_observer.hpp
 * Project: sdmmonitor
 * Created Date: Wednesday 8th May 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 21st May 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct to create observer which outputs condensation rate monitored from SDM microphysical
 * process in each gridbox a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS_SDMMONITOR_CONDENSATION_OBSERVER_HPP_
#define LIBS_OBSERVERS_SDMMONITOR_CONDENSATION_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>

#include "../../kokkosaliases.hpp"
#include "../consttstep_observer.hpp"
#include "../observers.hpp"
#include "./do_sdmmonitor_obs.hpp"

/* struct satisfies SDMMonitor concept for use in do_sdmmonitor_obs to make observer */
struct MonitorCondensation {
  using datatype = float;
  Buffer<datatype>::mirrorviewd_buffer d_data;  // must match view type used by DoSDMMonitorObs

  MonitorCondensation() : d_data("condrate", 1) { Kokkos::deep_copy(d_data, 0.0); }

  /**
   * @brief Monitor condensation rate
   *
   * _Note:_ conversion of condensation rate from double precision (8 bytes double) to single
   * precision (4 bytes float) in output.
   *
   */
  KOKKOS_FUNCTION
  void monitor_microphysics() const;  // TODO(CB) monitor condensation properly
};

/**
 * @brief Constructs an observer which writes data monitoring condensation microphysics to an
 * array with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @return Constructed type satisfying observer concept.
 */
template <typename Store>
inline Observer auto CondensationObserver(const unsigned int interval, Dataset<Store> &dataset,
                                          const size_t maxchunk) {
  const auto xzarr_ptr = std::make_shared<XarrayZarrArray<Store, MonitorCondensation::datatype>>(
      dataset.template create_array<MonitorCondensation::datatype>("condrate", "TODO(CB)", "<f4",
                                                                   0.5, {maxchunk}, {"time"}));

  const auto do_obs = DoSDMMonitorObs<Store, MonitorCondensation, MonitorCondensation::datatype>(
      dataset, xzarr_ptr);
  return ConstTstepObserver(interval, do_obs);
}

#endif  //  LIBS_OBSERVERS_SDMMONITOR_CONDENSATION_OBSERVER_HPP_
