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
 * Last Modified: Thursday 11th April 2024
 * Modified By: CB
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
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./collect_data_for_dataset.hpp"
#include "./generic_collect_data.hpp"
#include "./write_to_dataset_observer.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr/dataset.hpp"

template <typename Store>
struct RaggedCount {
 private:
  std::shared_ptr<XarrayZarrArray<Store, uint32_t>> xzarr_ptr;  ///< pointer to raggedcount array

 public:
  RaggedCount(const Dataset<Store> &dataset, const size_t maxchunk)
      : xzarr_ptr(std::make_shared<XarrayZarrArray<Store, uint32_t>>(
            dataset.template create_raggedcount_array<uint32_t>(
                "raggedcount", "", "<u4", 1, {maxchunk}, {"time"}, "superdroplets"))) {}

  /* writes the total number of super-droplets in the domain "totnsupers" to the raggedcount array
  in the dataset. Note static conversion from architecture dependent, usually 16 byte unsigned
  integer (size_t = uint64_t), to 8 byte unsigned integer (uint32_t). */
  void write_to_array(const Dataset<Store> &dataset, const viewd_constsupers totsupers) const {
    const auto totnsupers = static_cast<uint32_t>(totsupers.extent(0));
    dataset.write_to_array(xzarr_ptr, totnsupers);
  }

  void write_arrayshape(const Dataset<Store> &dataset) const {
    dataset.write_arrayshape(xzarr_ptr);
  }
};

/**
 * @brief Constructs type sastifying the CollectDataForDataset concept for a given Store (using an
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
 * @param maxchunk The maximum chunk size (number of elements).
 * @param ngbxs The number of gridboxes.
 * @return CollectDataForDataset<Store> An instance satisfying the CollectDataForDataset concept for
 * collecting a 2-D floating point variable (e.g. a thermodynamic variable) from each gridbox.
 */

/**
 * @brief Constructs type sastifying the CollectDataForDataset concept for a given Store (using an
 * instance of the GenericCollectData class) which writes a superdroplet variable to a ragged Xarray
 * in a dataset.
 *
 * Function return type writes a superdroplet varaible "name" to a ragged Xarray for a data
 * type by collecting data according to the given FunctorFunc from within a
 * Kokkos::parallel_for loop over superdroplets with a range policy.
 *
 * @tparam Store The type of the dataset store.
 * @tparam T The type of the superdroplet variable data.
 * @tparam FunctorFunc The type of the functor function.
 * @param dataset The dataset to collect data from.
 * @param ffunc The functor function to apply to collect data.
 * @param name The name of the variable.
 * @param units The units of the variable.
 * @param dtype A string representing the data type of the variable.
 * @param scale_factor The scale factor of the variable.
 * @param maxchunk The maximum chunk size.
 * @return CollectDataForDataset<Store> An instance of CollectDataForDataset for collecting
 * superdroplet variable data.
 */
template <typename Store, typename T, typename FunctorFunc>
CollectDataForDataset<Store> auto CollectSuperdropVariable(
    const Dataset<Store> &dataset, const FunctorFunc ffunc, const std::string_view name,
    const std::string_view units, const std::string_view dtype, const double scale_factor,
    const size_t maxchunk) {
  const auto chunkshape = std::vector<size_t>{maxchunk};
  const auto dimnames = std::vector<std::string>{"superdroplets"};
  const auto sampledimname = std::string_view("superdroplets");
  const auto xzarr = dataset.template create_ragged_array<T>(name, units, dtype, scale_factor,
                                                             chunkshape, dimnames, sampledimname);

  return GenericCollectData(ffunc, xzarr, 0);
}

/* Operator is functor to perform copy of value of superdroplet's gridbox index "sdgbxindex" for
each superdroplet in totsupers view to d_data in parallel. Note conversion of

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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' sdgbxindex.
 */
struct SdgbxindexFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto sdgbxindex = static_cast<uint32_t>(totsupers(kk).get_sdgbxindex());
    d_data(kk) = sdgbxindex;
  }
};

/**
 * @brief Constructs a collector for each superdroplets' gridbox index data.
 *
 * This function constructs a type satisfying the CollectDataForDataset to collect each
 * superdroplet's gridbox index "sdgbxindex" as a 4 byte unsigned integer and write it to a ragged
 * array in a dataset.
 *
 * @tparam Store The type of the dataset store.
 * @param dataset The dataset to collect data from.
 * @param maxchunk The maximum chunk size (number of elements).
 * @return CollectDataForDataset<Store> An instance of CollectDataForDataset for collecting
 * sdgbxindex data.
 */
template <typename Store>
CollectDataForDataset<Store> auto CollectSdgbxindex(const Dataset<Store> &dataset,
                                                    const int maxchunk) {
  return CollectSuperdropVariable<Store, uint32_t, SdgbxindexFunc>(
      dataset, SdgbxindexFunc{}, "sdgbxindex", "", "<u4", 1, maxchunk);
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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' identity.
 */
struct SdIdFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<uint32_t>::mirrorviewd_buffer d_data) const {
    auto sdid = static_cast<uint32_t>(totsupers(kk).sdId.value);
    d_data(kk) = sdid;
  }
};

template <typename Store>
CollectDataForDataset<Store> auto CollectSdId(const Dataset<Store> &dataset, const int maxchunk) {
  return CollectSuperdropVariable<Store, uint32_t, SdIdFunc>(dataset, SdIdFunc{}, "sdId", "", "<u4",
                                                             1, maxchunk);
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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' multiplicity.
 */
struct XiFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<uint64_t>::mirrorviewd_buffer d_data) const {
    assert((totsupers(kk).get_xi() < LIMITVALUES::uint64_t_max) &&
           "superdroplet mulitiplicy too large to represent with 4 byte unsigned integer");
    auto xi = static_cast<uint64_t>(totsupers(kk).get_xi());
    d_data(kk) = xi;
  }
};

template <typename Store>
CollectDataForDataset<Store> auto CollectXi(const Dataset<Store> &dataset, const int maxchunk) {
  return CollectSuperdropVariable<Store, uint64_t, XiFunc>(dataset, XiFunc{}, "xi", "", "<u8", 1,
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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' radius.
 */
struct RadiusFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto radius = static_cast<float>(totsupers(kk).get_radius());
    d_data(kk) = radius;
  }
};

template <typename Store>
CollectDataForDataset<Store> auto CollectRadius(const Dataset<Store> &dataset, const int maxchunk) {
  return CollectSuperdropVariable<Store, float, RadiusFunc>(
      dataset, RadiusFunc{}, "radius", "micro-m", "<f4", dlc::R0 * 1e6, maxchunk);
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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' msol.
 */
struct MsolFunc {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto msol = static_cast<float>(totsupers(kk).get_msol());
    d_data(kk) = msol;
  }
};

template <typename Store>
CollectDataForDataset<Store> auto CollectMsol(const Dataset<Store> &dataset, const int maxchunk) {
  return CollectSuperdropVariable<Store, float, MsolFunc>(dataset, MsolFunc{}, "msol", "g", "<f4",
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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' coord3.
 */
struct Coord3Func {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto coord3 = static_cast<float>(totsupers(kk).get_coord3());
    d_data(kk) = coord3;
  }
};

template <typename Store>
CollectDataForDataset<Store> auto CollectCoord3(const Dataset<Store> &dataset, const int maxchunk) {
  return CollectSuperdropVariable<Store, float, Coord3Func>(dataset, Coord3Func{}, "coord3", "m",
                                                            "<f4", dlc::COORD0, maxchunk);
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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' coord1.
 */
struct Coord1Func {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto coord1 = static_cast<float>(totsupers(kk).get_coord1());
    d_data(kk) = coord1;
  }
};

template <typename Store>
CollectDataForDataset<Store> auto CollectCoord1(const Dataset<Store> &dataset, const int maxchunk) {
  return CollectSuperdropVariable<Store, float, Coord1Func>(dataset, Coord1Func{}, "coord1", "m",
                                                            "<f4", dlc::COORD0, maxchunk);
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
 * @param totsupers The view of superdroplets on device.
 * @param d_data The mirror view buffer for the superdroplets' coord2.
 */
struct Coord2Func {
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t kk, viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
                  Buffer<float>::mirrorviewd_buffer d_data) const {
    auto coord2 = static_cast<float>(totsupers(kk).get_coord2());
    d_data(kk) = coord2;
  }
};

template <typename Store>
CollectDataForDataset<Store> auto CollectCoord2(const Dataset<Store> &dataset, const int maxchunk) {
  return CollectSuperdropVariable<Store, float, Coord2Func>(dataset, Coord2Func{}, "coord2", "m",
                                                            "<f4", dlc::COORD0, maxchunk);
}

/**
 * @brief Constructs an observer which writes superdroplet variables (e.g. their attributes) from
 * each superdroplet at start of each observation timestep to a ragged arrays with a constant
 * observation timestep "interval".
 *
 * @tparam Store Type of store for dataset.
 * @param interval Observation timestep.
 * @param dataset Dataset to write time data to.
 * @param maxchunk Maximum number of elements in a chunk (1-D vector size).
 * @param collect_data Object satisfying CollectDataForDataset for given Store to write superdroplet
 * data, such as their attributes, to ragged arrays.
 * @return Observer An observer instance for writing thermodynamic variables from each gridbox.
 */
template <typename Store>
inline Observer auto SuperdropsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                        const int maxchunk,
                                        CollectDataForDataset<Store> auto collect_data) {
  const CollectRaggedCount<Store> auto ragged_count = RaggedCount(dataset, maxchunk);
  return WriteToDatasetObserver(interval, dataset, collect_data, ragged_count);
}

#endif  // LIBS_OBSERVERS_SUPERDROPS_OBSERVER_HPP_
