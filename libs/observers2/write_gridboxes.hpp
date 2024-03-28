/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: write_gridboxes.hpp
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
 * Template for an observer which writes data collected from Gridboxes at the start of
 * each timestep in parallel to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_
#define LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/**
 * @brief Concept for GridboxDataWriter is all types that have functions for creating a functor to
 * collect data from gridboxes and then writing that data to an array in a dataset.
 *
 * @tparam GDW The type that satisfies the GridboxDataWriter concept.
 */
template <typename GDW, typename Store>
concept GridboxDataWriter = requires(GDW gdw, Dataset<Store> &ds, const viewd_constgbx d_gbxs) {
  { gdw.get_functor(d_gbxs) };
  { gdw.write_to_array(ds) } -> std::same_as<void>;
  { gdw.write_arrayshape(ds) } -> std::same_as<void>;
};

/**
 * @brief Combined gridbox data writer struct combines two gridbox data writers into one.
 *
 * @tparam GbxWriter1 The type of the first gridbox data writer.
 * @tparam GbxWriter2 The type of the second gridbox data writer.
 */
template <typename Store, GridboxDataWriter<Store> GbxWriter1, GridboxDataWriter<Store> GbxWriter2>
struct CombinedGridboxDataWriter {
 private:
  GbxWriter1 a; /**< The first instance of type of GridboxDataWriter. */
  GbxWriter2 b; /**< The second instance of type of GridboxDataWriter. */

 public:
  struct Functor {
    GbxWriter1::Functor a_functor;
    GbxWriter2::Functor b_functor;

    explicit Functor(const GbxWriter1 a, const GbxWriter2 b, const viewd_constgbx d_gbxs)
        : a_functor(a.get_functor(d_gbxs)), b_functor(b.get_functor(d_gbxs)) {}

    // Functor operator to perform copy of each element in parallel
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const {
      a_functor(ii);
      b_functor(ii);
    }
  };

  /**
   * @brief Constructs a CombinedGridboxDataWriter object.
   *
   * @param a The first gridbox data writer.
   * @param b The second gridbox data writer.
   */
  CombinedGridboxDataWriter(const GbxWriter1 a, const GbxWriter2 b) : a(a), b(b) {}

  Functor get_functor(const viewd_constgbx d_gbxs) const { return Functor(a, b, d_gbxs); }

  void write_to_array(Dataset<Store> &dataset) const {
    a.write_to_array(dataset);
    b.write_to_array(dataset);
  }

  void write_arrayshape(Dataset<Store> &dataset) const {
    a.write_arrayshape(dataset);
    b.write_arrayshape(dataset);
  }
};

/**
 * @brief Operator for combining two gridbox data writers which write to a FSStore.
 *
 * This operator combines two gridbox data writers into one using the
 * CombinedGridboxDataWriter struct with Store = FSStore.
 *
 * @param a The first gridbox data writer.
 * @param b The second gridbox data writer.
 * @return The combined gridbox data writer.
 */
template <typename Store>
struct CombineGDW {
  template <GridboxDataWriter<Store> GbxWriter1, GridboxDataWriter<Store> GbxWriter2>
  auto operator()(const GbxWriter1 a, const GbxWriter2 b) const {
    return CombinedGridboxDataWriter<Store, GbxWriter1, GbxWriter2>(a, b);
  }
};

// struct satifying GridboxDataWriter and does nothing
template <typename Store>
struct NullGbxWriter {
 public:
  struct Functor {
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const {}
  };

  Functor get_functor(const viewd_constgbx d_gbxs) const { return Functor{}; }

  void write_to_array(Dataset<Store> &dataset) const {}

  void write_arrayshape(Dataset<Store> &dataset) const {}
};

// template GridboxDataWriter to write one variable from each gridbox to an array in a dataset
template <typename Store, typename T, typename FunctorFunc>
class OneVarGbxWriter {
 private:
  FunctorFunc ffunc;
  using viewh_data = Buffer<T>::viewh_buffer;              // type of view for h_data
  using mirrorviewd_data = Buffer<T>::mirrorviewd_buffer;  // mirror view type for d_data
  std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr;    // pointer to array in dataset
  viewh_data h_data;        // view on host for value of 1 variable from every gridbox
  mirrorviewd_data d_data;  // mirror view of h_data on device

 public:
  template <typename FunctorFunc>
  struct Functor {
    FunctorFunc ffunc;
    viewd_constgbx d_gbxs;    // view of gridboxes
    mirrorviewd_data d_data;  // mirror view on device for value of 1 variable from every gridbox

    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, mirrorviewd_data d_data)
        : ffunc(ffunc), d_gbxs(d_gbxs), d_data(d_data) {}

    // Functor operator to perform copy of 1 variable in each gridbox to d_data in parallel
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const { ffunc(ii) }
  };

  // Constructor to initialize views and pointer to array in dataset
  OneVarGbxWriter(Dataset<Store> &dataset, FunctorFunc ffunc,
                  std::make_shared<XarrayZarrArray<Store, T>> xzarr_ptr, const size_t ngbxs)
      : ffunc(ffunc),
        xzarr_ptr(xzarr_ptr),
        h_data("h_data", ngbxs),
        d_data(Kokkos::create_mirror_view(ExecSpace(), h_data)) {}

  // return functor for getting pressure from each gridbox in parallel
  Functor get_functor(const viewd_constgbx d_gbxs) const {
    assert((d_gbxs.extent(0) == d_data.extent(0)) &&
           "d_data view must be size of the number of gridboxes");
    return Functor(ffunc, d_gbxs, d_data);
  }

  // copy data from device view directly to host and then write to array in dataset
  void write_to_array(Dataset<Store> &dataset) const {
    Kokkos::deep_copy(h_data, d_data);
    dataset.write_to_array(xzarr_ptr, h_data);
  }

  // call function to write shape of array according to dataset
  void write_arrayshape(Dataset<Store> &dataset) const { dataset.write_arrayshape(xzarr_ptr); }
};

/* template class for observing variables from each gridbox in parallel
and then writing them to their repspective arrays in a dataset */
template <typename Store, GridboxDataWriter<Store> GbxWriter>
class WriteGridboxes {
 private:
  Dataset<Store> &dataset;  ///< dataset to write data to
  GbxWriter writer;  ///< object collects data from gridboxes and writes it to arrays in the dataset

  // use functor from writer to collect data from gridboxes in parallel
  void collect_data_from_gridboxes(const viewd_constgbx d_gbxs) const {
    const size_t ngbxs(d_gbxs.extent(0));
    auto functor = writer.get_functor(d_gbxs);
    Kokkos::parallel_for("collect_gbxs_data", Kokkos::RangePolicy<ExecSpace>(0, ngbxs), functor);
  }

  // collect data from gridboxes and then write it to arrays in the dataset
  void at_start_step(const viewd_constgbx d_gbxs) const {
    collect_data_from_gridboxes(d_gbxs);
    writer.write_to_array(dataset);
    // dataset.set_dimension({"time", time+1}); // TODO(CB) do this with coord observer
  }

 public:
  WriteGridboxes(Dataset<Store> &dataset, GbxWriter writer) : dataset(dataset), writer(writer) {}

  ~WriteGridboxes() { writer.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes Gridboxes to Dataset observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    at_start_step(d_gbxs);
  }
};

#endif  // LIBS_OBSERVERS2_WRITE_GRIDBOXES_HPP_
