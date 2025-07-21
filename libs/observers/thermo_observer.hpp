/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: thermo_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to write variables related to Gridboxes' state at the start of
 * a constant interval timestep to arrays in a dataset
 */

#ifndef LIBS_OBSERVERS_THERMO_OBSERVER_HPP_
#define LIBS_OBSERVERS_THERMO_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "observers/collect_data_for_dataset.hpp"
#include "observers/generic_collect_data.hpp"
#include "observers/write_to_dataset_observer.hpp"

/**
 * @brief Constructs type sastifying the CollectDataForDataset concept for a given Dataset (using an
 * instance of the GenericCollectData class) which writes a thermodynamic variable to an Xarray in a
 * dataset.
 *
 * Function return type writes a thermodyanmic varaible "name" to an Xarray as a 4-byte floating
 * point type by collecting data according to the given FunctorFunc from within a
 * Kokkos::parallel_for loop over gridboxes with a range policy.
 *
 * @param dataset The dataset to write the variable to.
 * @param ffunc The functor function to collect the variable from within a parallel range policy
 * over gridboxes.
 * @param name The name of the new coordinate.
 * @param units The units of the coordinate.
 * @param scale_factor The scale factor of the coordinate data.
 * @param maxchunk The maximum chunk size (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return CollectDataForDataset<Dataset> An instance satisfying the CollectDataForDataset concept
 * for collecting a 2-D floating point variable (e.g. a thermodynamic variable) from each gridbox.
 */
template <typename Dataset, typename FunctorFunc>
CollectDataForDataset<Dataset> auto CollectThermoVariable(
    const Dataset &dataset, const FunctorFunc ffunc, const std::string_view name,
    const std::string_view units, const double scale_factor, const size_t maxchunk,
    const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  const auto xzarr =
      dataset.template create_array<float>(name, units, scale_factor, chunkshape, dimnames);
  return GenericCollectData(ffunc, xzarr, ngbxs);
}

/**
 * @brief Functor operator to perform a copy of the pressure from the state of each gridbox to
 * d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of press from double (8 bytes) to single precision float (4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the pressure in each gridbox.
 */
struct PressFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto press = static_cast<float>(d_gbxs(ii).state.press);
    d_data(ii) = press;
  }
};

/**
 * @brief Functor operator to perform a copy of the temperature from the state of each gridbox to
 * d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of temp from double (8 bytes) to single precision float (4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the temperature in each gridbox.
 */
struct TempFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto temp = static_cast<float>(d_gbxs(ii).state.temp);
    d_data(ii) = temp;
  }
};

/**
 * @brief Functor operator to perform a copy of the vapour mass mixing ratio "qvap" from the state
 * of each gridbox to d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of qvap from double (8 bytes) to single precision float (4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for qvap from each gridbox.
 */
struct QvapFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto qvap = static_cast<float>(d_gbxs(ii).state.qvap);
    d_data(ii) = qvap;
  }
};

/**
 * @brief Functor operator to perform a copy of the liquid mass mixing ratio "qcond" from the state
 * of each gridbox to d_data within Kokkos::parallel_for loop over gridboxes with range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of qcond from double (8 bytes) to single precision float (4 bytes).
 *
 * @param ii The index of the gridbox.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for qcond from each gridbox.
 */
struct QcondFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto qcond = static_cast<float>(d_gbxs(ii).state.qcond);
    d_data(ii) = qcond;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept for collecting multiple
 * thermodyanmcis variables from each gridbox and writing them to a dataset.
 *
 * This function combines CollectDataForDataset types for many thermodynamic variables from each
 * gridbox (e.g. press, temp, qvap, qcond, etc.) using instances of the GenericCollectData class.
 *
 * @tparam Dataset The type of dataset.
 * @param dataset The dataset to write the wind velocity components to.
 * @param maxchunk The maximum chunk size (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * thermodynamics from the state of each gridbox.
 */
template <typename Dataset>
inline CollectDataForDataset<Dataset> auto CollectThermo(const Dataset &dataset,
                                                         const size_t maxchunk,
                                                         const size_t ngbxs) {
  const CollectDataForDataset<Dataset> auto press = CollectThermoVariable<Dataset, PressFunc>(
      dataset, PressFunc{}, "press", "hPa", dlc::P0 / 100, maxchunk, ngbxs);

  const CollectDataForDataset<Dataset> auto temp = CollectThermoVariable<Dataset, TempFunc>(
      dataset, TempFunc{}, "temp", "K", dlc::TEMP0, maxchunk, ngbxs);

  const CollectDataForDataset<Dataset> auto qvap = CollectThermoVariable<Dataset, QvapFunc>(
      dataset, QvapFunc{}, "qvap", "g/Kg", 1000.0, maxchunk, ngbxs);

  const CollectDataForDataset<Dataset> auto qcond = CollectThermoVariable<Dataset, QcondFunc>(
      dataset, QcondFunc{}, "qcond", "g/Kg", 1000.0, maxchunk, ngbxs);

  const CollectDataForDataset<Dataset> auto a = CombinedCollectDataForDataset(press, temp);
  const CollectDataForDataset<Dataset> auto b = CombinedCollectDataForDataset(qvap, qcond);
  return CombinedCollectDataForDataset(a, b);
}

/**
 * @brief Constructs an observer which writes thermodyanmcis from each gridbox (e.g. press, temp,
 * qvap, etc.) at start of each observation timestep to a arrays with a constant
 * observation timestep "interval".
 *
 * @tparam Dataset Type of dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Observer An observer instance for writing thermodynamic variables from each gridbox.
 */
template <typename Dataset>
inline Observer auto ThermoObserver(const unsigned int interval, const Dataset &dataset,
                                    const size_t maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Dataset> auto thermo = CollectThermo(dataset, maxchunk, ngbxs);
  return WriteToDatasetObserver(interval, dataset, thermo);
}

#endif  // LIBS_OBSERVERS_THERMO_OBSERVER_HPP_
