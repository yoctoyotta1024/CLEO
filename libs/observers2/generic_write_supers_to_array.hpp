/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: generic_write_supers_to_array.hpp
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
 * Concept and structs to write data collected from all the superdrops in the domain ("totsupers")
 * in parallel to a ragged array in a dataset.
 */

#ifndef LIBS_OBSERVERS2_GENERIC_WRITE_SUPERS_TO_ARRAY_HPP_
#define LIBS_OBSERVERS2_GENERIC_WRITE_SUPERS_TO_ARRAY_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>
#include <iostream>
#include <memory>

#include "../kokkosaliases.hpp"
#include "./write_gridbox_to_array.hpp"
#include "./xarray_for_supers_data.hpp"  // TODO(CB)
#include "superdrops/superdrop.hpp"
#include "zarr2/dataset.hpp"

/* template WriteGridboxToArray to write one variable from all the superdrops in the domain,
"totsupers", to a ragged array in a dataset */
template <typename Store, typename T, typename FunctorFunc>
class GenericWriteSupersToXarray {
 private:
  std::shared_ptr<XarrayForSupersData<Store, T>> xzarr_ptr;
  FunctorFunc ffunc;

 public:
  struct Functor {
    using mirrorviewd_data = XarrayForSupersData<Store, T>::mirrorviewd_data;
    FunctorFunc ffunc;
    viewd_constsupers totsupers;  // view on device of the superdroplets
    mirrorviewd_data d_data;      // mirror view for data on device

    Functor(FunctorFunc ffunc, const viewd_constsupers totsupers, mirrorviewd_data d_data)
        : ffunc(ffunc), totsupers(totsupers), d_data(d_data) {}

    /* Functor operator to perform copy of 1 variable in each gridbox to d_data in parallel when
    using Kokkos range policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const size_t kk) const { ffunc(kk, totsupers, d_data); }
  };

  /* Constructor to initialize views and pointer to array in dataset */
  GenericWriteSupersToXarray(const Dataset<Store> &dataset, const std::string_view name,
                             const std::string_view units, const std::string_view dtype,
                             const double scale_factor, const size_t maxchunk, FunctorFunc ffunc)
      : xzarr_ptr(std::make_shared<XarrayForSupersData<Store, T>>(dataset, name, units, dtype,
                                                                  scale_factor, maxchunk)),
        ffunc(ffunc) {}

  /* return functor for getting 1 variable from every superdroplets in each gridbox in parallel */
  Functor get_functor(const viewd_constsupers totsupers) const {
    Kokkos::realloc(xzarr_ptr->h_data, totsupers.extent(0));  // TODO(CB) move into XarraySupersData
    Kokkos::realloc(xzarr_ptr->d_data, totsupers.extent(0));
    assert((totsupers.extent(0) == d_data.extent(0)) &&
           "d_data view must be size of the total number of superdroplets");
    return Functor(ffunc, totsupers, xzarr_ptr->d_data);
  }

  // copy data from device view directly to host and then write to array in dataset
  void write_to_array(const Dataset<Store> &dataset) const { xzarr_ptr->write_to_array(dataset); }

  // call function to write shape of array according to dataset
  void write_arrayshape(const Dataset<Store> &dataset) const {
    xzarr_ptr->write_arrayshape(dataset);
  }
};

#endif  // LIBS_OBSERVERS2_GENERIC_WRITE_SUPERS_TO_ARRAY_HPP_
