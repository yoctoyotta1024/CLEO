/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: generic_write_gridbox_to_array.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 3rd April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Concept and structs to write data collected from each Gridbox in parallel
 * to an array in a dataset
 */

#ifndef LIBS_OBSERVERS2_GENERIC_WRITE_GRIDBOX_TO_ARRAY_HPP_
#define LIBS_OBSERVERS2_GENERIC_WRITE_GRIDBOX_TO_ARRAY_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "./write_gridbox_to_array.hpp"
#include "./xarray_for_gridbox_data.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr2/dataset.hpp"

/* template WriteGridboxToArray to write one variable from each gridbox to an array in a dataset */
template <typename Store, typename T, typename FunctorFunc>
class GenericWriteGridboxToXarray {
 private:
  std::shared_ptr<XarrayForGridboxData<Store, T>> xzarr_ptr;
  FunctorFunc ffunc;

 public:
  struct Functor {
    using mirrorviewd_data = XarrayForGridboxData<Store, T>::mirrorviewd_data;
    FunctorFunc ffunc;
    viewd_constgbx d_gbxs;    // view of gridboxes on device
    mirrorviewd_data d_data;  // mirror view for data on device

    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, mirrorviewd_data d_data)
        : ffunc(ffunc), d_gbxs(d_gbxs), d_data(d_data) {}

    /* Functor operator to perform copy of 1 variable in each gridbox to d_data in parallel when
    using Kokkos range policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t ii) const { ffunc(ii, d_gbxs, d_data); }

    /* Functor operator to perform copy of 1 variable in each gridbox to d_data in parallel
    when using Kokkos range policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const { ffunc(team_member, d_gbxs, d_data); }
  };

  /* Constructor to initialize views and pointer to array in dataset */
  GenericWriteGridboxToXarray(const Dataset<Store> &dataset, const std::string_view name,
                              const std::string_view units, const std::string_view dtype,
                              const double scale_factor, const size_t maxchunk, const size_t ngbxs,
                              FunctorFunc ffunc)
      : xzarr_ptr(std::make_shared<XarrayForGridboxData<Store, T>>(dataset, name, units, dtype,
                                                                   scale_factor, maxchunk, ngbxs)),
        ffunc(ffunc) {}

  /* return functor for getting 1 variable from every gridbox in parallel */
  Functor get_functor(const viewd_constgbx d_gbxs) const {
    assert((d_gbxs.extent(0) == xzarr_ptr->d_data.extent(0)) &&
           "d_data view must be size of the number of gridboxes");
    return Functor(ffunc, d_gbxs, xzarr_ptr->d_data);
  }

  /* copy data from device view directly to host and then write to array in dataset */
  void write_to_array(const Dataset<Store> &dataset) const { xzarr_ptr->write_to_array(dataset); }

  /* call function to write shape of array according to dataset */
  void write_arrayshape(const Dataset<Store> &dataset) const {
    xzarr_ptr->write_arrayshape(dataset);
  }
};

#endif  // LIBS_OBSERVERS2_GENERIC_WRITE_GRIDBOX_TO_ARRAY_HPP_
