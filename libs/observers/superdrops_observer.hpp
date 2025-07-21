/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: superdrops_observer.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to write variables related to superdroplet attributes at the start of
 * a constant interval timestep to ragged arrays in a dataset
 */

#ifndef LIBS_OBSERVERS_SUPERDROPS_OBSERVER_HPP_
#define LIBS_OBSERVERS_SUPERDROPS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "observers/collect_data_for_dataset.hpp"
#include "observers/generic_collect_data.hpp"
#include "observers/write_to_dataset_observer.hpp"
#include "superdrops/superdrop.hpp"

/**
 * @brief A struct for collecting ragged count data.
 *
 * This struct is responsible for collecting ragged count data, which represents the count of the
 * number of super-droplets written during a write of a ragged array of superdroplet data.
 *
 * @tparam Dataset Type of dataset.
 * @tparam Store Type of store for dataset.
 */
template <typename Dataset, typename Store>
struct RaggedCount {
 private:
  std::shared_ptr<XarrayZarrArray<Store, uint32_t>> xzarr_ptr; /**< pointer to raggedcount Xarray */

 public:
  /**
   * @brief Constructs a RaggedCount object.
   *
   * Constructs a RaggedCount object with the specified dataset and maximum chunk size.
   *
   * @param dataset The dataset to collect data from.
   * @param store The store the dataset writes to.
   * @param maxchunk The maximum chunk size (number of elements).
   */
  RaggedCount(const Dataset &dataset, Store &store, const size_t maxchunk)
      : xzarr_ptr(std::make_shared<XarrayZarrArray<Store, uint32_t>>(
            dataset.template create_raggedcount_array<uint32_t>("raggedcount", "", 1, {maxchunk},
                                                                {"time"}, "superdroplets"))) {}

  /**
   * @brief Writes the total number of super-droplets to the ragged count array in the dataset.
   *
   * _Note:_ static conversion from architecture dependent, usually 16 byte unsigned integer (size_t
   * = uint64_t), to 8 byte unsigned integer (uint32_t).
   *
   * @param dataset The dataset to write data to.
   * @param d_supers The view of total super-droplets.
   */
  void write_to_array(const Dataset &dataset, const subviewd_constsupers d_supers) const {
    const auto totnsupers = static_cast<uint32_t>(d_supers.extent(0));
    dataset.write_to_array(xzarr_ptr, totnsupers);
  }

  /**
   * @brief Writes the shape of the ragged count array to the dataset.
   *
   * This function writes the shape of the ragged count array to the dataset.
   *
   * @param dataset The dataset to write data to.
   */
  void write_arrayshape(const Dataset &dataset) const { dataset.write_arrayshape(xzarr_ptr); }
};

/**
 * @brief Constructs type sastifying the CollectDataForDataset concept for a given Dataset (using an
 * instance of the GenericCollectData class) which writes a superdroplet variable to a ragged Xarray
 * in a dataset.
 *
 * Function return type writes a superdroplet varaible "name" to a ragged Xarray for a data
 * type by collecting data according to the given FunctorFunc from within a
 * Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * @tparam Dataset The type of the dataset.
 * @tparam T The type of the superdroplet variable data.
 * @tparam FunctorFunc The type of the functor function.
 * @param dataset The dataset to collect data from.
 * @param ffunc The functor function to apply to collect data.
 * @param name The name of the variable.
 * @param units The units of the variable.
 * @param scale_factor The scale factor of the variable.
 * @param maxchunk The maximum chunk size.
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * superdroplet variable data.
 */
template <typename Dataset, typename T, typename FunctorFunc>
CollectDataForDataset<Dataset> auto CollectSuperdropVariable(
    const Dataset &dataset, const FunctorFunc ffunc, const std::string_view name,
    const std::string_view units, const double scale_factor, const size_t maxchunk) {
  const auto chunkshape = std::vector<size_t>{maxchunk};
  const auto dimnames = std::vector<std::string>{"superdroplets"};
  const auto sampledimname = std::string_view("superdroplets");
  const auto xzarr = dataset.template create_ragged_array<T>(name, units, scale_factor, chunkshape,
                                                             dimnames, sampledimname);

  return GenericCollectData(ffunc, xzarr, 0);
}

/**
 * @brief Functor operator to perform a copy of the sdgbxindex of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of sdgbxindex from unsigned long int (8 bytes) to unsigned int
 * (uint32_t = 4 bytes)
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' sdgbxindex.
 */
struct SdgbxindexFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto sdgbxindex = static_cast<uint32_t>(d_supers(kk).get_sdgbxindex());
    d_data(kk) = sdgbxindex;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept for each superdroplets'
 * gridbox index data.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's gridbox index "sdgbxindex" as a 4 byte unsigned integer and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * sdgbxindex data.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectSdgbxindex(const Dataset &dataset,
                                                      const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, uint32_t, SdgbxindexFunc>(dataset, SdgbxindexFunc{},
                                                                     "sdgbxindex", "", 1, maxchunk);
}

/**
 * @brief Functor operator to perform a copy of the identity of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of sdid from unsigned long int (8 bytes) to unsigned int
 * (uint32_t = 4 bytes)
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' identity.
 */
struct SdIdFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto sdid = static_cast<uint32_t>(d_supers(kk).sdId.get_value());
    d_data(kk) = sdid;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset conceptfor each superdroplets'
 * identity.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's identity "sdId" as a 4 byte unsigned integer and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * sdId data.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectSdId(const Dataset &dataset, const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, uint32_t, SdIdFunc>(dataset, SdIdFunc{}, "sdId", "", 1,
                                                               maxchunk);
}

/**
 * @brief Functor operator to perform a copy of the multiplicity of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of xi from size_t (architecture dependent usually 8 bytes) to 8 byte, long
 * unsigned int (unit64_t).
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' multiplicity.
 */
struct XiFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<uint64_t>::mirrorviewd_buffer d_data) const {
    assert((d_supers(kk).get_xi() < LIMITVALUES::uint64_t_max) &&
           "superdroplet mulitiplicy too large to represent with 4 byte unsigned integer");
    auto xi = static_cast<uint64_t>(d_supers(kk).get_xi());
    d_data(kk) = xi;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept for each superdroplets'
 * multiplicity.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's multiplicity "xi" as a 8 byte unsigned integer and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting xi
 * data.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectXi(const Dataset &dataset, const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, uint64_t, XiFunc>(dataset, XiFunc{}, "xi", "", 1,
                                                             maxchunk);
}

/**
 * @brief Functor operator to perform a copy of the radius of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of radius from double (8 bytes) to float (4 bytes).
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' radius.
 */
struct RadiusFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto radius = static_cast<float>(d_supers(kk).get_radius());
    d_data(kk) = radius;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept or each superdroplets'
 * radius.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's radius as a 4 byte floating point value and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * radius.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectRadius(const Dataset &dataset, const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, float, RadiusFunc>(dataset, RadiusFunc{}, "radius",
                                                              "micro-m", dlc::R0 * 1e6, maxchunk);
}

/**
 * @brief Functor operator to perform a copy of the solute mass of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of msol from double (8 bytes) to float (4 bytes).
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' msol.
 */
struct MsolFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto msol = static_cast<float>(d_supers(kk).get_msol());
    d_data(kk) = msol;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept or each superdroplets'
 * solute mass.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's solute mass "msol" as a 4 byte floating point value and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting msol.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectMsol(const Dataset &dataset, const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, float, MsolFunc>(dataset, MsolFunc{}, "msol", "g",
                                                            dlc::MASS0grams, maxchunk);
}

/**
 * @brief Functor operator to perform a copy of the coord3 of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of coord3 from double (8 bytes) to float (4 bytes).
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' coord3.
 */
struct Coord3Func {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto coord3 = static_cast<float>(d_supers(kk).get_coord3());
    d_data(kk) = coord3;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept or each superdroplets'
 * coord3.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's coord3 as a 4 byte floating point value and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * coord3.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectCoord3(const Dataset &dataset, const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, float, Coord3Func>(dataset, Coord3Func{}, "coord3", "m",
                                                              dlc::COORD0, maxchunk);
}

/**
 * @brief Functor operator to perform a copy of the coord1 of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of coord1 from double (8 bytes) to float (4 bytes).
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' coord1.
 */
struct Coord1Func {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto coord1 = static_cast<float>(d_supers(kk).get_coord1());
    d_data(kk) = coord1;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept or each superdroplets'
 * coord1.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's coord1 as a 4 byte floating point value and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * coord1.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectCoord1(const Dataset &dataset, const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, float, Coord1Func>(dataset, Coord1Func{}, "coord1", "m",
                                                              dlc::COORD0, maxchunk);
}

/**
 * @brief Functor operator to perform a copy of the coord2 of each superdroplet to
 * d_data within Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * Signature of operator such that type can be used by GenericCollectData struct for FunctorFunc.
 *
 * _Note:_ Conversion of coord2 from double (8 bytes) to float (4 bytes).
 *
 * @param kk The index of the superdrop.
 * @param d_gbxs The view of gridboxes on device.
 * @param d_supers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' coord2.
 */
struct Coord2Func {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const subviewd_constsupers d_supers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto coord2 = static_cast<float>(d_supers(kk).get_coord2());
    d_data(kk) = coord2;
  }
};

/**
 * @brief Constructs a type satisyfing the CollectDataForDataset concept or each superdroplets'
 * coord2.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's coord2 as a 4 byte floating point value and write it to a ragged
 * array in a dataset.
 *
 * @tparam Dataset The type of the dataset.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Dataset> An instance of CollectDataForDataset for collecting
 * coord2.
 */
template <typename Dataset>
CollectDataForDataset<Dataset> auto CollectCoord2(const Dataset &dataset, const size_t maxchunk) {
  return CollectSuperdropVariable<Dataset, float, Coord2Func>(dataset, Coord2Func{}, "coord2", "m",
                                                              dlc::COORD0, maxchunk);
}

/**
 * @brief Constructs an observer which writes superdroplet variables (e.g. their attributes) from
 * each superdroplet at start of each observation timestep to a ragged arrays with a constant
 * observation timestep "interval".
 *
 * @tparam Dataset Type of dataset.
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param store The store the dataset writes to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param collect_data Object satisfying CollectDataForDataset for given Store to write superdroplet
 * data, such as their attributes, to ragged arrays.
 * @return Observer An observer instance for writing thermodynamic variables from each gridbox.
 */
template <typename Dataset, typename Store>
inline Observer auto SuperdropsObserver(const unsigned int interval, const Dataset &dataset,
                                        Store &store, const size_t maxchunk,
                                        CollectDataForDataset<Dataset> auto collect_data) {
  const CollectRaggedCount<Dataset> auto ragged_count = RaggedCount(dataset, store, maxchunk);
  return WriteToDatasetObserver(interval, dataset, collect_data, ragged_count);
}

#endif  // LIBS_OBSERVERS_SUPERDROPS_OBSERVER_HPP_
