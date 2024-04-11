/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: massmoments_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 11th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output the mass moments of the droplet size distribution in each gridbox
 * to individual arrays in a dataset a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS_MASSMOMENTS_OBSERVER_HPP_
#define LIBS_OBSERVERS_MASSMOMENTS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <concepts>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./collect_data_for_dataset.hpp"
#include "./generic_collect_data.hpp"
#include "./observers.hpp"
#include "./write_to_dataset_observer.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr/dataset.hpp"

/**
 * @brief Functor to perform calculation of 0th, 1st, and 2nd moments of the (real)
 * droplet mass distribution in each gridbox.
 *
 * This operator is a functor to perform the calculation of the 0th, 1st, and 2nd moments
 * of the droplet mass distribution in each gridbox within a Kokkos::parallel_for range policy
 * loop over superdroplets. The calculation is equivalent to 0th, 3rd, and 6th moments of the
 * droplet radius distribution for each gridbox.
 *
 * _Note:_ the conversion from 8 to 4-byte precision for all three moments.
 */
struct MassMomentsFunc {
  /**
   * @brief Functor operator to perform calculation of 0th, 1st, and 2nd moments of the (real)
   * droplet mass distribution in each gridbox.
   *
   * This operator is a functor to perform the calculation of the 0th, 1st, and 2nd moments
   * of the droplet mass distribution in each gridbox (i.e. 0th, 3rd, and 6th moments of the
   * droplet radius distribution) within a Kokkos::parallel_reduce range policy
   * loop over superdroplets.
   *
   * Kokkos::parallel_reduce([...]) is equivalent in serial to sum over result of:
   * for (size_t kk(0); kk < supers.extent(0); ++kk){[...]}.
   *
   * _Note:_ conversion from 8 to 4-byte precision for all mass moments: mom0 from size_t
   * (architecture dependent usually long unsigned int = 8 bytes) to 8 byte unsigned integer, and
   * mom1 and mom2 from double (8 bytes) to float (4 bytes).
   * @param team_member The Kokkos team member.
   * @param d_gbxs The view of gridboxes on device.
   * @param d_mom0 The mirror view buffer for the 0th mass moment.
   * @param d_mom1 The mirror view buffer for the 1st mass moment.
   * @param d_mom2 The mirror view buffer for the 2nd mass moment.
   */
  KOKKOS_FUNCTION
  void operator()(const TeamMember &team_member, const viewd_constgbx d_gbxs,
                  Buffer<uint64_t>::mirrorviewd_buffer d_mom0,
                  Buffer<float>::mirrorviewd_buffer d_mom1,
                  Buffer<float>::mirrorviewd_buffer d_mom2) const;
};

/**
 * @brief Functor to perform calculation of 0th, 1st, and 2nd moments of the (real)
 * rain-droplet mass distribution in each gridbox.
 *
 * This operator is a functor to perform the calculation of the 0th, 1st, and 2nd moments
 * of the rain-droplet mass distribution in each gridbox within a Kokkos::parallel_for range
 * policy loop over superdroplets. The calculation is equivalent to 0th, 3rd, and 6th moments of the
 * rain-droplet radius distribution for each gridbox.
 *
 * A raindroplet is a droplet with a radius >= rlim = 40microns.
 *
 * _Note:_ the conversion from 8 to 4-byte precision for all three moments.
 */
struct RaindropsMassMomentsFunc {
  /**
   * @brief Functor operator to perform calculation of 0th, 1st, and 2nd moments of the (real)
   * droplet mass distribution in each gridbox.
   *
   * This operator is a functor to perform the calculation of the 0th, 1st, and 2nd moments
   * of the droplet mass distribution in each gridbox (i.e. 0th, 3rd, and 6th moments of the
   * droplet radius distribution) within a Kokkos::parallel_for range policy
   * loop over superdroplets.
   *
   * A raindroplet is a droplet with a radius >= rlim = 40microns.
   *
   * Kokkos::parallel_reduce([...]) is equivalent in serial to sum over result of:
   * for (size_t kk(0); kk < supers.extent(0); ++kk){[...]}.
   *
   * _Note:_ conversion from 8 to 4-byte precision for all mass moments: mom0 from size_t
   * (architecture dependent usually long unsigned int = 8 bytes) to 8 byte unsigned integer, and
   * mom1 and mom2 from double (8 bytes) to float (4 bytes).
   *
   * @param team_member The Kokkos team member.
   * @param d_gbxs The view of gridboxes on device.
   * @param d_mom0 The mirror view buffer for the 0th mass moment.
   * @param d_mom1 The mirror view buffer for the 1st mass moment.
   * @param d_mom2 The mirror view buffer for the 2nd mass moment.
   */
  KOKKOS_FUNCTION
  void operator()(const TeamMember &team_member, const viewd_constgbx d_gbxs,
                  Buffer<uint64_t>::mirrorviewd_buffer d_mom0,
                  Buffer<float>::mirrorviewd_buffer d_mom1,
                  Buffer<float>::mirrorviewd_buffer d_mom2) const;
};

/**
 * @brief A struct satisfying CollectDataForDataset for collecting the 0th, 1st, and 2nd moments of
 * the (rain)droplet mass distribution in each gridbox.
 *
 * This struct is similar to the GenericCollectData struct but specialized with xarrays and a
 * functor that stores three variables (the 0th, 1st, and 2nd mass moments) and with a functor which
 * can act inside a Kokkos team policy rather than a range policy.
 */
template <typename Store, typename FunctorFunc>
struct CollectMassMoments {
 private:
  FunctorFunc ffunc;
  std::shared_ptr<XarrayAndViews<Store, uint64_t>> mom0_ptr; /**< Xarray for 0th mass moment */
  std::shared_ptr<XarrayAndViews<Store, float>> mom1_ptr;    /**< Xarray for 1st mass moment */
  std::shared_ptr<XarrayAndViews<Store, float>> mom2_ptr;    /**< Xarray for 2nd mass moment */

  /**
   * @brief Deep copies data from a device view to the host and then writes the host data to an
   * array in the dataset.
   *
   * @tparam T The data type.
   * @param ptr A pointer to a struct with the data views to be copied and the Xarray to write to.
   * @param dataset The dataset for the Xarray.
   */
  template <typename T>
  void write_one_array(std::shared_ptr<XarrayAndViews<Store, T>> ptr,
                       const Dataset<Store> &dataset) const {
    Kokkos::deep_copy(ptr->h_data, ptr->d_data);
    dataset.write_to_array(ptr->xzarr, ptr->h_data);
  }

  /**
   * @brief Calls a function to write the shape of the array according to the dataset.
   *
   * @tparam T The data type.
   * @param ptr A pointer to struct with the Xarray whose shape must be written.
   * @param dataset The dataset to write array shape to.
   */
  template <typename T>
  void write_one_arrayshape(std::shared_ptr<XarrayAndViews<Store, T>> ptr,
                            const Dataset<Store> &dataset) const {
    dataset.write_arrayshape(ptr->xzarr);
  }

 public:
  /**
   * @brief Generic functor to collect all three mass moments from within a Kokkos::parallel_for
   * loop over gridboxes with a team policy.
   */
  struct Functor {
    using mirrorviewd_mom0 = XarrayAndViews<Store, uint64_t>::mirrorviewd_data;
    using mirrorviewd_mom1 = XarrayAndViews<Store, float>::mirrorviewd_data;
    using mirrorviewd_mom2 = XarrayAndViews<Store, float>::mirrorviewd_data;
    FunctorFunc ffunc; /**< Functor to calculate mass moments within parallel team policy loop */
    viewd_constgbx d_gbxs;   /**< View of gridboxes on device */
    mirrorviewd_mom0 d_mom0; /**< Mirror view on device for 0th mass moment */
    mirrorviewd_mom1 d_mom1; /**< Mirror view on device for 1st mass moment */
    mirrorviewd_mom2 d_mom2; /**< Mirror view on device for 2nd mass moment */

    /**
     * @brief Constructs a Functor object.
     *
     * @param ffunc The functor to calculate mass moments within a parallel team policy loop.
     * @param d_gbxs The view of gridboxes on the device.
     * @param d_mom0 The mirror view for the 0th mass moment on the device.
     * @param d_mom1 The mirror view for the 1st mass moment on the device.
     * @param d_mom2 The mirror view for the 2nd mass moment on the device.
     */
    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, mirrorviewd_mom0 d_mom0,
            mirrorviewd_mom1 d_mom1, mirrorviewd_mom2 d_mom2)
        : ffunc(ffunc), d_gbxs(d_gbxs), d_mom0(d_mom0), d_mom1(d_mom1), d_mom2(d_mom2) {}

    /**
     * @brief Adapter to call functor to perform calculation of massmoments in each gridbox
     * from within a Kokkos::parallel_for loop with a team policy.
     */
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {
      ffunc(team_member, d_gbxs, d_mom0, d_mom1, d_mom2);
    }
  };

  /**
   * @brief Constructor to initialize CollectMassMoments.
   *
   * @param ffunc Functor function-like object for calculting mass moments within parallel loop.
   * @param xzarr_mom0 Xarray for the 0th mass moment.
   * @param xzarr_mom1 Xarray for the 1st mass moment.
   * @param xzarr_mom2 Xarray for the 2nd mass moment.
   * @param dataview_size Size of the data view used to collect data from within the functor
   * function call.
   */
  CollectMassMoments(const FunctorFunc ffunc, const XarrayZarrArray<Store, uint64_t> xzarr_mom0,
                     const XarrayZarrArray<Store, float> xzarr_mom1,
                     const XarrayZarrArray<Store, float> xzarr_mom2, const size_t dataview_size)
      : ffunc(ffunc),
        mom0_ptr(std::make_shared<XarrayAndViews<Store, uint64_t>>(xzarr_mom0, dataview_size)),
        mom1_ptr(std::make_shared<XarrayAndViews<Store, float>>(xzarr_mom1, dataview_size)),
        mom2_ptr(std::make_shared<XarrayAndViews<Store, float>>(xzarr_mom2, dataview_size)) {}

  /**
   * @brief Returns a functor for getting 0th, 1st, and 2nd mass moments from every gridbox
   * within a Kokkos::parallel_for loop over gridboxes with a team policy.
   *
   * @param d_gbxs View of gridboxes on device.
   * @param totsupers View of superdroplets on device.
   * @return Functor for collecting mass moments.
   */
  Functor get_functor(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    assert(((mom0_ptr->d_data.extent(0) == d_gbxs.extent(0)) &&
            (mom1_ptr->d_data.extent(0) == d_gbxs.extent(0)) &&
            (mom2_ptr->d_data.extent(0) == d_gbxs.extent(0))) &&
           "d_data views for mass moments should be size of the number of gridboxes");
    return Functor(ffunc, d_gbxs, mom0_ptr->d_data, mom1_ptr->d_data, mom2_ptr->d_data);
  }

  /**
   * @brief Writes all three mass moments to arrays in the dataset.
   *
   * @param dataset The dataset to write data to.
   */
  void write_to_arrays(const Dataset<Store> &dataset) const {
    write_one_array(mom0_ptr, dataset);
    write_one_array(mom1_ptr, dataset);
    write_one_array(mom2_ptr, dataset);
  }

  /**
   * @brief Writes the shape of all three arrays to the dataset.
   *
   * @param dataset The dataset to write data to.
   */
  void write_arrayshapes(const Dataset<Store> &dataset) const {
    write_one_arrayshape(mom0_ptr, dataset);
    write_one_arrayshape(mom1_ptr, dataset);
    write_one_arrayshape(mom2_ptr, dataset);
  }

  /**
   * @brief Null function to satisfy CollectDataForDataset concept.
   *
   * @param dataset The dataset to write data to.
   */
  void write_to_ragged_arrays(const Dataset<Store> &dataset) const {}

  /**
   * @brief Null function to satisfy CollectDataForDataset concept.
   *
   * @param dataset The dataset to write data to.
   */
  void write_ragged_arrayshapes(const Dataset<Store> &dataset) const {}

  /**
   * @brief Null function to satisfy CollectDataForDataset concept.
   *
   * @param sz The size for reallocation of a view.
   */
  void reallocate_views(const size_t sz) const {}
};

/**
 * @brief Creates an XarrayZarrArray for storing the mass moments of each gridbox in a dataset.
 *
 * @tparam Store The type of data store in the dataset.
 * @tparam T The type of the mass moment data to store in the XarrayZarrArray.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param name The name of the Xarray.
 * @param units The units of the mass moment data.
 * @param dtype A string representing the data type of the Xarray.
 * @param scale_factor The scale factor for the data.
 * @param maxchunk The maximum chunk size (number of elements) for the Xarray.
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, T> The created XarrayZarrArray.
 */
template <typename Store, typename T>
XarrayZarrArray<Store, T> create_massmoment_xarray(const Dataset<Store> &dataset,
                                                   const std::string_view name,
                                                   const std::string_view units,
                                                   const std::string_view dtype,
                                                   const double scale_factor, const size_t maxchunk,
                                                   const size_t ngbxs) {
  const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
  const auto dimnames = std::vector<std::string>{"time", "gbxindex"};
  return dataset.template create_array<T>(name, units, dtype, scale_factor, chunkshape, dimnames);
}

/**
 * @brief Creates an XarrayZarrArray for storing the 0th mass moment in a dataset.
 *
 * Calls create_massmoment_xarray for data that is represented by 8 byte unsigned integers with
 * no units and is called "name" - e.g. the 0th mass moment of a droplet distribution.
 *
 * @tparam Store The type of data store in the dataset.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param name The name of the (0th mass moment) Xarray.
 * @param maxchunk The maximum chunk size for the Xarray (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, uint64_t> The created XarrayZarrArray (for the 0th mass moment).
 */
template <typename Store>
XarrayZarrArray<Store, uint64_t> create_massmom0_xarray(const Dataset<Store> &dataset,
                                                        const std::string_view name,
                                                        const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("");
  return create_massmoment_xarray<Store, uint64_t>(dataset, name, units, "<u8", 1, maxchunk, ngbxs);
}

/**
 * @brief Creates an XarrayZarrArray for storing the 1st mass moment in a dataset.
 *
 * Calls create_massmoment_xarray for data that is represented by 4 byte float with
 * units "g" and is called "name" - e.g. the 1st mass moment of a droplet distribution.
 *
 * @tparam Store The type of data store in the dataset.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param name The name of the (1st mass moment) Xarray.
 * @param maxchunk The maximum chunk size for the Xarray (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, uint64_t> The created XarrayZarrArray (for the 1st mass moment).
 */
template <typename Store>
XarrayZarrArray<Store, float> create_massmom1_xarray(const Dataset<Store> &dataset,
                                                     const std::string_view name,
                                                     const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("g");
  constexpr auto scale_factor = dlc::MASS0grams;
  return create_massmoment_xarray<Store, float>(dataset, name, units, "<f4", scale_factor, maxchunk,
                                                ngbxs);
}

/**
 * @brief Creates an XarrayZarrArray for storing the 2nd mass moment in a dataset.
 *
 * Calls create_massmoment_xarray for data that is represented by 4 byte float with
 * units "g^2" and is called "name" - e.g. the 2nd mass moment of a droplet distribution.
 *
 * @tparam Store The type of data store in the dataset.
 * @param dataset The dataset where the XarrayZarrArray will be created.
 * @param name The name of the (2nd mass moment) Xarray.
 * @param maxchunk The maximum chunk size for the Xarray (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return XarrayZarrArray<Store, uint64_t> The created XarrayZarrArray (for the 2nd mass moment).
 */
template <typename Store>
XarrayZarrArray<Store, float> create_massmom2_xarray(const Dataset<Store> &dataset,
                                                     const std::string_view name,
                                                     const size_t maxchunk, const size_t ngbxs) {
  const auto units = std::string_view("g^2");
  constexpr auto scale_factor = dlc::MASS0grams * dlc::MASS0grams;
  return create_massmoment_xarray<Store, float>(dataset, name, units, "<f4", scale_factor, maxchunk,
                                                ngbxs);
}

/**
 * @brief Constructs an observer which writes mass moments of droplet distribution at start of
 * each observation timestep to an array with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Constructed type for writing mass moments of droplet distribution satisfying the
 * observer concept.
 */
template <typename Store>
inline Observer auto MassMomentsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                         const int maxchunk, const size_t ngbxs) {
  const auto xzarr_mom0 = create_massmom0_xarray(dataset, "massmom0", maxchunk, ngbxs);
  const auto xzarr_mom1 = create_massmom1_xarray(dataset, "massmom1", maxchunk, ngbxs);
  const auto xzarr_mom2 = create_massmom2_xarray(dataset, "massmom2", maxchunk, ngbxs);

  const auto ffunc = MassMomentsFunc{};

  const CollectDataForDataset<Store> auto massmoments =
      CollectMassMoments(ffunc, xzarr_mom0, xzarr_mom1, xzarr_mom2, ngbxs);
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesTeamPolicyFunc{}, dataset, massmoments);
  return WriteToDatasetObserver(interval, parallel_write);
}

/**
 * @brief Constructs an observer which writes mass moments of rain-droplet distribution at start
 * of each observation timestep to an array with a constant observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param ngbxs The number of gridboxes.
 * @return Constructed type for writing mass moments of rain-droplet distribution satisfying the
 * observer concept.
 */
template <typename Store>
inline Observer auto MassMomentsRaindropsObserver(const unsigned int interval,
                                                  const Dataset<Store> &dataset, const int maxchunk,
                                                  const size_t ngbxs) {
  const auto xzarr_mom0 = create_massmom0_xarray(dataset, "massmom0_raindrops", maxchunk, ngbxs);
  const auto xzarr_mom1 = create_massmom1_xarray(dataset, "massmom1_raindrops", maxchunk, ngbxs);
  const auto xzarr_mom2 = create_massmom2_xarray(dataset, "massmom2_raindrops", maxchunk, ngbxs);

  const auto ffunc = RaindropsMassMomentsFunc{};

  const CollectDataForDataset<Store> auto massmoments_raindrops =
      CollectMassMoments(ffunc, xzarr_mom0, xzarr_mom1, xzarr_mom2, ngbxs);
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesTeamPolicyFunc{}, dataset, massmoments_raindrops);
  return WriteToDatasetObserver(interval, parallel_write);
}

#endif  // LIBS_OBSERVERS_MASSMOMENTS_OBSERVER_HPP_
