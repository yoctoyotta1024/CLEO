/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: create_massmoments_arrays.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Helper functions to create arrays for mass moments of the droplet size distribution in each
 * gridbox at each timestep interval.
 */

#ifndef LIBS_OBSERVERS_CREATE_MASSMOMENTS_ARRAYS_HPP_
#define LIBS_OBSERVERS_CREATE_MASSMOMENTS_ARRAYS_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <concepts>
#include <cstdint>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "observers/collect_data_for_dataset.hpp"
#include "observers/generic_collect_data.hpp"
#include "observers/observers.hpp"
#include "observers/write_to_dataset_observer.hpp"
#include "superdrops/superdrop.hpp"

/**
 * @brief Creates an XarrayZarrArray for storing the mass moments of each gridbox in a dataset.
 *
 * @tparam Dataset The type of dataset.
 * @tparam Store The type of data store in the dataset.
 * @tparam T The type of the mass moment data to store in the XarrayZarrArray.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param store The store the dataset writes to.
 * @param name The name of the Xarray.
 * @param units The units of the mass moment data.
 * @param scale_factor The scale factor for the data.
 * @param maxchunk The maximum chunk size (number of elements) for the Xarray.
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, T> The created XarrayZarrArray.
 */
template <typename Dataset, typename Store, typename T>
XarrayZarrArray<Store, T> create_massmoment_xarray(const Dataset &dataset, Store &store,
                                                   const std::string_view name,
                                                   const std::string_view units,
                                                   const double scale_factor, const size_t maxchunk,
                                                   const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  return dataset.template create_array<T>(name, units, scale_factor, chunkshape, dimnames);
}

/**
 * @brief Creates an XarrayZarrArray for storing the 0th mass moment in a dataset.
 *
 * Calls create_massmoment_xarray for data that is represented by 8 byte unsigned integers with
 * no units and is called "name" - e.g. the 0th mass moment of a droplet distribution.
 *
 * @tparam Dataset The type of dataset.
 * @tparam Store The type of data store in the dataset.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param store The store the dataset writes to.
 * @param name The name of the (0th mass moment) Xarray.
 * @param maxchunk The maximum chunk size for the Xarray (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, uint64_t> The created XarrayZarrArray (for the 0th mass moment).
 */
template <typename Dataset, typename Store>
XarrayZarrArray<Store, uint64_t> create_massmom0_xarray(const Dataset &dataset, Store &store,
                                                        const std::string_view name,
                                                        const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("");
  return create_massmoment_xarray<Dataset, Store, uint64_t>(dataset, store, name, units, 1,
                                                            maxchunk, ngbxs);
}

/**
 * @brief Creates an XarrayZarrArray for storing the 1st mass moment in a dataset.
 *
 * Calls create_massmoment_xarray for data that is represented by 4 byte float with
 * units "g" and is called "name" - e.g. the 1st mass moment of a droplet distribution.
 *
 * @tparam Dataset The type of dataset.
 * @tparam Store The type of data store in the dataset.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param store The store the dataset writes to.
 * @param name The name of the (1st mass moment) Xarray.
 * @param maxchunk The maximum chunk size for the Xarray (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, uint64_t> The created XarrayZarrArray (for the 1st mass moment).
 */
template <typename Dataset, typename Store>
XarrayZarrArray<Store, float> create_massmom1_xarray(const Dataset &dataset, Store &store,
                                                     const std::string_view name,
                                                     const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("g");
  constexpr auto scale_factor = dlc::MASS0grams;
  return create_massmoment_xarray<Dataset, Store, float>(dataset, store, name, units, scale_factor,
                                                         maxchunk, ngbxs);
}

/**
 * @brief Creates an XarrayZarrArray for storing the 2nd mass moment in a dataset.
 *
 * Calls create_massmoment_xarray for data that is represented by 4 byte float with
 * units "g^2" and is called "name" - e.g. the 2nd mass moment of a droplet distribution.
 *
 * @tparam Dataset The type of dataset.
 * @tparam Store The type of data store in the dataset.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param store The store the dataset writes to.
 * @param name The name of the (2nd mass moment) Xarray.
 * @param maxchunk The maximum chunk size for the Xarray (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, uint64_t> The created XarrayZarrArray (for the 2nd mass moment).
 */
template <typename Dataset, typename Store>
XarrayZarrArray<Store, float> create_massmom2_xarray(const Dataset &dataset, Store &store,
                                                     const std::string_view name,
                                                     const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("g^2");
  constexpr auto scale_factor = dlc::MASS0grams * dlc::MASS0grams;
  return create_massmoment_xarray<Dataset, Store, float>(dataset, store, name, units, scale_factor,
                                                         maxchunk, ngbxs);
}

#endif  // LIBS_OBSERVERS_CREATE_MASSMOMENTS_ARRAYS_HPP_
