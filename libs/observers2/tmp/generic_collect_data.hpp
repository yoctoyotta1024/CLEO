/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: generic_collect_data.hpp
 * Project: tmp
 * Created Date: Thursday 4th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 4th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * generic struct satisyfing the CollectDataForDataset concept to collect data a
 * variable(s) and write it to xarray(s) in a datatset.
 */

#ifndef LIBS_OBSERVERS2_TMP_GENERIC_COLLECT_DATA_HPP_
#define LIBS_OBSERVERS2_TMP_GENERIC_COLLECT_DATA_HPP_

#include <Kokkos_Core.hpp>
#include <memory>

#include "../kokkosaliases.hpp"
#include "zarr2/buffer.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"

/* generic struct satisyfing the CollectDataForDataset concept to collect data for a
 * variable and write it to an xarray in a datatset. */
template <typename Store, typename T, typename FunctorFunc>
class GenericCollectData {
 private:
  using viewh_data = Buffer<T>::viewh_buffer;              // type of view for h_data
  using mirrorviewd_data = Buffer<T>::mirrorviewd_buffer;  // mirror view type for d_data
  FunctorFunc ffunc;
  std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr;
  viewh_data h_data;        // view on host for value of 1 variable from every superdrop
  mirrorviewd_data d_data;  // mirror view of h_data on device

 public:
  struct Functor {
    FunctorFunc ffunc;
    viewd_constgbx d_gbxs;        // view of gridboxes on device
    viewd_constsupers totsupers;  // view of superdroplets on device
    mirrorviewd_data d_data;      // mirror view for data to collect on device

    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, const viewd_constsupers totsupers,
            mirrorviewd_data d_data)
        : ffunc(ffunc), d_gbxs(d_gbxs), totsupers(totsupers), d_data(d_data) {}

    /* Functor operator to perform copy of 1 variable in gridboxes and/or superdroplets
    to d_data from within a parallel loop using a Kokkos range policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t nn) const { ffunc(nn, d_gbxs, totsupers, d_data); }

    /* Functor operator to perform copy of 1 variable in gridboxes and/or superdroplets
    to d_data from within a parallel loop using a Kokkos team policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {
      ffunc(team_member, d_gbxs, totsupers, d_data);
    }
  };

  /* Constructor to initialize GenericCollectData given functor function-like object, shared pointer
  to an xarray in a dataset and the size of the data view used to collect data
  from within the functor function call. */
  GenericCollectData(const FunctorFunc ffunc,
                     const std::shared_ptr<XarrayZarrArray<Store, T>> xzarr_ptr,
                     const size_t dataview_size)
      : ffunc(ffunc),
        xzarr_ptr(xzarr_ptr),
        h_data("h_data", dataview_size),
        d_data(Kokkos::create_mirror_view(ExecSpace(), h_data)) {}

  /* return functor for getting 1 variable from every gridbox in parallel */
  Functor get_functor(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    assert(((d_gbxs.extent(0) == d_data.extent(0)) || (totsupers.extent(0) == d_data.extent(0))) &&
           "d_data view should be size of the number of gridboxes or superdroplets");
    return Functor(ffunc, d_gbxs, totsupers, d_data);
  }

  void reallocate_dataviews(const size_t size) {
    // TODO(CB) fix forseen issues over const-ness
    Kokkos::realloc(h_data, size);
    Kokkos::realloc(d_data, size);
  }

  /* copy data from device view directly to host and then write to array in dataset */
  void write_to_arrays(const Dataset<Store> &dataset) const {
    Kokkos::deep_copy(h_data, d_data);
    dataset.write_to_array(xzarr_ptr, h_data);
  }

  /* call function to write shape of array according to dataset */
  void write_arrayshapes(const Dataset<Store> &dataset) const {
    dataset.write_arrayshape(xzarr_ptr);
  }
};

#endif  // LIBS_OBSERVERS2_TMP_GENERIC_COLLECT_DATA_HPP_
