/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: collect_data_for_dataset.hpp
 * Project: observers
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Concept of CollectDataForDataset and monoidal structs which can be used within called to
 * ParallelWriteData operator() to collect data within parallel loops and write it to arrays in a
 * dataset
 */

#ifndef LIBS_OBSERVERS_COLLECT_DATA_FOR_DATASET_HPP_
#define LIBS_OBSERVERS_COLLECT_DATA_FOR_DATASET_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../kokkosaliases.hpp"

/**
 * @brief Concept for CollectDataForDataset is all types that have functions for creating a functor
 * to collect data from a gridbox and/or superdroplets (to use in a Kokkos parallel loop) and then
 * write the data to arrays in a dataset.
 *
 * @tparam CDD The type that satisfies the CollectDataForDataset concept.
 */
template <typename CDD, typename Dataset>
concept CollectDataForDataset = requires(CDD cdd, const Dataset &ds, const viewd_constgbx d_gbxs,
                                         const subviewd_constsupers d_supers, const size_t sz) {
  { cdd.get_functor(d_gbxs, d_supers) };
  { cdd.reallocate_views(sz) } -> std::same_as<void>;
  { cdd.template write_to_arrays<Dataset>(ds) } -> std::same_as<void>;
  { cdd.template write_to_ragged_arrays<Dataset>(ds) } -> std::same_as<void>;
  { cdd.template write_arrayshapes<Dataset>(ds) } -> std::same_as<void>;
  { cdd.template write_ragged_arrayshapes<Dataset>(ds) } -> std::same_as<void>;
};

/**
 * @brief struct is a new type satirfying CollectDataForDataset formed from the combination of
 * two structs that also satisfy the CollectDataForDataset concept (given the same Dataset type).
 * This new type does the actions of the original structs in sequence.
 *
 * Note the contrains of the CollectDataForDataset<Dataset> are not applied directly on the incoming
 * two types (CollectData1 and CollectData2) that are combined. They should each nevertheless
 * obey the contraints of the CollectDataForDataset<Dataset> concept given the same Dataset type.
 *
 * @tparam CollectData1 The type of the first CollectDataForDataset<Dataset>.
 * @tparam CollectData2 The type of the second CollectDataForDataset<Dataset>.
 */
template <typename CollectData1, typename CollectData2>
struct CombinedCollectDataForDataset {
 private:
  CollectData1 a; /**< The first instance of type of CollectDataForDataset. */
  CollectData2 b; /**< The second instance of type of CollectDataForDataset. */

 public:
  struct Functor {
    CollectData1::Functor a_functor;
    CollectData2::Functor b_functor;

    explicit Functor(const CollectData1 a, const CollectData2 b, const viewd_constgbx d_gbxs,
                     const subviewd_constsupers d_supers)
        : a_functor(a.get_functor(d_gbxs, d_supers)), b_functor(b.get_functor(d_gbxs, d_supers)) {}

    /* Functor operator to perform copy of each element in parallel in Kokkos Range Policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t nn) const {
      a_functor(nn);
      b_functor(nn);
    }

    /* Functor operator to perform copy of each element in parallel in Kokkos Team Policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {
      a_functor(team_member);
      b_functor(team_member);
    }
  };

  /**
   * @brief Constructs a CombinedCollectDataForDataset object.
   *
   * @param a The first type of CollectDataForDataset instance.
   * @param b The second type of CollectDataForDataset instance.
   */
  CombinedCollectDataForDataset(const CollectData1 a, const CollectData2 b) : a(a), b(b) {}

  Functor get_functor(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    return Functor(a, b, d_gbxs, d_supers);
  }

  template <typename Dataset>
  void write_to_arrays(const Dataset &dataset) const {
    a.write_to_arrays(dataset);
    b.write_to_arrays(dataset);
  }

  template <typename Dataset>
  void write_to_ragged_arrays(const Dataset &dataset) const {
    a.write_to_ragged_arrays(dataset);
    b.write_to_ragged_arrays(dataset);
  }

  template <typename Dataset>
  void write_arrayshapes(const Dataset &dataset) const {
    a.write_arrayshapes(dataset);
    b.write_arrayshapes(dataset);
  }

  template <typename Dataset>
  void write_ragged_arrayshapes(const Dataset &dataset) const {
    a.write_ragged_arrayshapes(dataset);
    b.write_ragged_arrayshapes(dataset);
  }

  void reallocate_views(const size_t sz) const {
    a.reallocate_views(sz);
    b.reallocate_views(sz);
  }
};

/* struct satifying CollectDataForDataset and does nothing */
struct NullCollectDataForDataset {
 public:
  struct Functor {
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t nn) const {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {}
  };

  Functor get_functor(const viewd_constgbx d_gbxs, const subviewd_constsupers d_supers) const {
    return Functor{};
  }

  template <typename Dataset>
  void write_to_arrays(const Dataset &dataset) const {}

  template <typename Dataset>
  void write_to_ragged_arrays(const Dataset &dataset) const {}

  template <typename Dataset>
  void write_arrayshapes(const Dataset &dataset) const {}

  template <typename Dataset>
  void write_ragged_arrayshapes(const Dataset &dataset) const {}

  void reallocate_views(const size_t sz) const {}
};

#endif  // LIBS_OBSERVERS_COLLECT_DATA_FOR_DATASET_HPP_
