/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: massmoments_observer.hpp
 * Project: observers2
 * Created Date: Wednesday 24th January 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 5th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output the mass moments of the droplet size distribution in each gridbox
 * to individual arrays in a dataset a constant interval at the start of each timestep.
 */

#ifndef LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <concepts>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./collect_data_for_dataset.hpp"
#include "./generic_collect_data.hpp"
#include "./observers.hpp"
#include "./write_to_dataset_observer.hpp"
#include "gridboxes/gridbox.hpp"
#include "superdrops/superdrop.hpp"
#include "zarr2/dataset.hpp"

/* struct satifying CollectDataForDataset for collecting the 0th, 1st and 2nd moments of the
 * (rain)droplet mass distribution in each gridbox (ie. 0th, 3rd and 6th moments of the radius
 * distribution). struct similar to GenericCollectData but specialised with xarrays and a functor
 * that stores 3 variables (0th, 1st and 2nd mass moments) and with functor which can act inside a
 * kokkos team policy not range policy (see signature of functor's operator() function) */
template <typename Store>
struct CollectMassMoments {
 private:
  FunctorFunc ffunc;
  std::shared_ptr<XarrayAndViews<Store, uint32_t>> mom0_ptr;  // xarray struct for 0th mass moment
  std::shared_ptr<XarrayAndViews<Store, float>> mom1_ptr;     // xarray struct for 1st mass moment
  std::shared_ptr<XarrayAndViews<Store, float>> mom2_ptr;     // xarray struct for 2nd mass moment

  /* copy data from device view directly to host and then write to array in dataset */
  template <typename T>
  void write_one_array(std::shared_ptr<XarrayAndViews<Store, T>> ptr,
                       const Dataset<Store> &dataset) const {
    Kokkos::deep_copy(ptr->h_data, ptr->d_data);
    dataset.write_to_array(ptr->xzarr, ptr->h_data);
  }

  /* call function to write shape of array according to dataset */
  template <typename T>
  void write_one_arrayshape(std::shared_ptr<XarrayAndViews<Store, T>> ptr,
                            const Dataset<Store> &dataset) const {
    dataset.write_arrayshape(ptr->xzarr);
  }

 public:
  /* functor to collect 3 variables from within a parallel team policy */
  struct Functor {
    using mirrorviewd_mom0 = XarrayAndViews<Store, uint32_t>::mirrorviewd_data;
    using mirrorviewd_mom1 = XarrayAndViews<Store, float>::mirrorviewd_data;
    using mirrorviewd_mom2 = XarrayAndViews<Store, float>::mirrorviewd_data;
    FunctorFunc ffunc;        // functor to calculate mass moments within parallel team policy loop
    viewd_constgbx d_gbxs;    // view of gridboxes on device
    mirrorviewd_mom0 d_mom0;  // mirror view 0th mass moment on device
    mirrorviewd_mom1 d_mom1;  // mirror view 1st mass moment on device
    mirrorviewd_mom2 d_mom2;  // mirror view 2nd mass moment on device

    Functor(FunctorFunc ffunc, const viewd_constgbx d_gbxs, mirrorviewd_mom0 d_mom0,
            mirrorviewd_mom0 d_mom1, mirrorviewd_mom0 d_mom2)
        : ffunc(ffunc), d_gbxs(d_gbxs), d_mom0(d_mom0), d_mom1(d_mom1), d_mom2(d_mom2) {}

    /* Functor operator to perform calculation of massmoments in each gridbox
    and then copy to d_data from within a parallel loop using a Kokkos team policy */
    KOKKOS_INLINE_FUNCTION
    void operator()(const TeamMember &team_member) const {
      ffunc(team_member, d_gbxs, d_mom0, d_mom1, d_mom2);
    }
  };

  /* Constructor to initialize CollectMassMoments given functor function-like object,
  the xarrays for the 0th, 1st and 2nd mass moments in the dataset and the size of the data view
  used to collect data from within the functor function call. */
  CollectMassMoments(const FunctorFunc ffunc, const XarrayZarrArray<Store, uint32_t> xzarr_mom0,
                     const XarrayZarrArray<Store, float> xzarr_mom1,
                     const XarrayZarrArray<Store, float> xzarr_mom2, const size_t dataview_size)
      : ffunc(ffunc),
        mom0_ptr(std::make_shared<XarrayAndViews<Store, uint32_t>>(xzarr_mom0, dataview_size)),
        mom1_ptr(std::make_shared<XarrayAndViews<Store, float>>(xzarr_mom1, dataview_size)),
        mom2_ptr(std::make_shared<XarrayAndViews<Store, float>>(xzarr_mom2, dataview_size)) {}

  /* return functor for getting 0th, 1st and 2nd mass moments from every gridbox in parallel team
   * policy */
  Functor get_functor(const viewd_constgbx d_gbxs, const viewd_constsupers totsupers) const {
    assert(((mom0_ptr->d_data.extent(0) == d_gbxs.extent(0)) &&
            (mom1_ptr->d_data.extent(0) == d_gbxs.extent(0)) &&
            (mom2_ptr->d_data.extent(0) == d_gbxs.extent(0))) &&
           "d_data views for mass moments should be size of the number of gridboxes");
    return Functor(ffunc, d_gbxs, mom0_ptr->d_data, mom1_ptr->d_data, mom2_ptr->d_data);
  }

  void write_to_arrays(const Dataset<Store> &dataset) const {
    write_one_array(mom0_ptr, dataset);
    write_one_array(mom1_ptr, dataset);
    write_one_array(mom2_ptr, dataset);
  }

  void write_arrayshapes(const Dataset<Store> &dataset) const {
    write_one_arrayshape(mom0_ptr, dataset);
    write_one_arrayshape(mom1_ptr, dataset);
    write_one_arrayshape(mom2_ptr, dataset);
  }

  void write_to_ragged_arrays(const Dataset<Store> &dataset) const {}

  void write_ragged_arrayshapes(const Dataset<Store> &dataset) const {}

  void reallocate_views(const size_t sz) const {}
};

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the WriteToDatasetObserver class */
template <typename Store>
inline Observer auto MassMomentsObserver(const unsigned int interval, const Dataset<Store> &dataset,
                                         const int maxchunk, const size_t ngbxs) {
  const CollectDataForDataset<Store> auto massmoments =
      CollectMassMoments(ffunc, xzarr_mom0, xzarr_mom1, xzarr_mom2, ngbxs);  // TODO(CB) WIP
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesTeamPolicyFunc{}, dataset, massmoments);
  return WriteToDatasetObserver(interval, parallel_write);
}

/* constructs observer which writes mass moments of raindroplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the WritetoDatasetObserver class */
template <typename Store>
inline Observer auto MassMomentsRaindropsObserver(const unsigned int interval,
                                                  const Dataset<Store> &dataset, const int maxchunk,
                                                  const size_t ngbxs) {
  const CollectDataForDataset<Store> auto massmoments_raindrops =
      CollectMassMoments(ffunc, xzarr_mom0, xzarr_mom1, xzarr_mom2, ngbxs);  // TODO(CB) WIP
  const auto parallel_write =
      ParallelWriteGridboxes(ParallelGridboxesTeamPolicyFunc{}, dataset, massmoments_raindrops);
  return WriteToDatasetObserver(interval, parallel_write);
}

#endif  // LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_
