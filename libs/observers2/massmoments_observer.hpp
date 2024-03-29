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
 * Last Modified: Friday 29th March 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Observer to output variables related to Gridboxes' state at the start of
 * each timestep to individual arrays in a dataset
 */

#ifndef LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_
#define LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_

#include <Kokkos_Core.hpp>
#include <array>
#include <concepts>
#include <iostream>
#include <memory>
#include <string_view>

#include "../cleoconstants.hpp"
#include "./observers.hpp"
#include "./write_gridbox_to_array.hpp"
#include "zarr2/buffer.hpp"
#include "zarr2/dataset.hpp"
#include "zarr2/xarray_zarr_array.hpp"
#include "zarr2/zarr_array.hpp"

/* calculated 0th, 1st and 2nd moments of the (real)
droplet mass distribution for one gridbox, i.e. 0th, 3rd and 6th
moments of the droplet radius distribution for one gridbox.
Kokkos::parallel_reduce([...]) is equivalent in serial to:
for (size_t kk(0); kk < supers.extent(0); ++kk){[...]} */
KOKKOS_FUNCTION
void calculate_massmoments(const TeamMember &team_member, const int ii,
                           const subviewd_constsupers supers,
                           Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                           Buffer<float>::mirrorviewd_buffer d_mom1,
                           Buffer<float>::mirrorviewd_buffer d_mom2);

/* calculated 0th, 1st and 2nd moments of the (real)
droplet mass distribution in each gridbox, i.e. 0th, 3rd and 6th
moments of the droplet radius distribution for each gridbox.
Calculation is done for all gridboxes in parallel.
Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}  */
void calculate_massmoments(const viewd_constgbx d_gbxs, Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                           Buffer<float>::mirrorviewd_buffer d_mom1,
                           Buffer<float>::mirrorviewd_buffer d_mom2);

/* calculated 0th, 1st and 2nd moments of the (real)
raindroplet mass distribution for one gridbox, i.e. 0th, 3rd and 6th
moments of the raindroplet radius distribution for one gridbox.
A raindrop is droplet with a radius >= rlim = 40microns.
Kokkos::parallel_reduce([...]) is equivalent in serial to:
for (size_t kk(0); kk < supers.extent(0); ++kk){[...]} */
KOKKOS_FUNCTION
void calculate_massmoments_raindrops(const TeamMember &team_member, const int ii,
                                     const subviewd_constsupers supers,
                                     Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                                     Buffer<float>::mirrorviewd_buffer d_mom1,
                                     Buffer<float>::mirrorviewd_buffer d_mom2);

/* calculated 0th, 1st and 2nd moments of the (real)
raindroplet mass distribution in each gridbox, i.e. 0th, 3rd and 6th
moments of the raindroplet radius distribution for each gridbox.
A raindrop is droplet with a radius >= rlim = 40microns.
Calculation is done for all gridboxes in parallel.
Kokkos::parallel_for([...]) is equivalent in serial to:
for (size_t ii(0); ii < d_gbxs.extent(0); ++ii){[...]}  */
void calculate_massmoments_raindrops(const viewd_constgbx d_gbxs,
                                     Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                                     Buffer<float>::mirrorviewd_buffer d_mom1,
                                     Buffer<float>::mirrorviewd_buffer d_mom2);

template <typename Store>
struct MassMomArrays {
  XarrayForGenericGbxWriter<Store, uint32_t> mom0_xzarr;  ///< 0th mass moment array
  XarrayForGenericGbxWriter<Store, float> mom1_xzarr;     ///< 1st mass moment array
  XarrayForGenericGbxWriter<Store, float> mom2_xzarr;     ///< 2nd mass moment array

  MassMomArrays(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs,
                std::array<std::string_view, 3> names)
      : mom0_xzarr(dataset, names.at(0), "", "<u4", 1, maxchunk, ngbxs),
        mom1_xzarr(dataset, names.at(1), "g", "<f4", dlc::MASS0grams, maxchunk, ngbxs),
        mom2_xzarr(dataset, names.at(2), "g^2", "<f4", dlc::MASS0grams * dlc::MASS0grams, maxchunk,
                   ngbxs) {}

  void write_arrayshape(Dataset<Store> &dataset) {
    mom0_xzarr.write_arrayshape(dataset);
    mom1_xzarr.write_arrayshape(dataset);
    mom2_xzarr.write_arrayshape(dataset);
  }

  void write_massmoments_to_arrays(Dataset<Store> &dataset) {
    mom0_xzarr.write_to_array(dataset);
    mom1_xzarr.write_to_array(dataset);
    mom2_xzarr.write_to_array(dataset);
  }

  XarrayForGenericGbxWriter<Store, uint32_t>::mirrorviewd_data get_d_mom0() {
    return mom0_xzarr.d_data;
  }

  XarrayForGenericGbxWriter<Store, float>::mirrorviewd_data get_d_mom1() {
    return mom1_xzarr.d_data;
  }

  XarrayForGenericGbxWriter<Store, float>::mirrorviewd_data get_d_mom2() {
    return mom2_xzarr.d_data;
  }
};

/* template class for observing 0th, 1st and 2nd mass moment (i.e. 0th, 3rd and 6th radius moment)
of droplet distribution in each gridbox to arrays in a dataset in a store */
template <typename Store, typename MassMomsCalc>
class DoMassMomsObs {
 private:
  Dataset<Store> &dataset;                           ///< dataset to write moments to
  std::shared_ptr<MassMomArrays<Store>> xzarrs_ptr;  ///< pointer to mass moment arrays in dataset
  MassMomsCalc calculate_massmoments;  ///< function like object to perform moment calculations

 public:
  DoMassMomsObs(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs,
                const std::array<std::string_view, 3> &names, MassMomsCalc calculate_massmoments)
      : dataset(dataset),
        xzarrs_ptr(std::make_shared<MassMomArrays<Store>>(dataset, maxchunk, ngbxs, names)),
        calculate_massmoments(calculate_massmoments) {}

  ~DoMassMomsObs() { xzarrs_ptr->write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes time observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    auto d_mom0 = xzarrs_ptr->get_d_mom0();
    auto d_mom1 = xzarrs_ptr->get_d_mom1();
    auto d_mom2 = xzarrs_ptr->get_d_mom2();

    calculate_massmoments(d_gbxs, d_mom0, d_mom1, d_mom2);
    xzarrs_ptr->write_massmoments_to_arrays(dataset);
  }
};

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto MassMomentsObserver(const unsigned int interval, Dataset<Store> &dataset,
                                         const int maxchunk, const size_t ngbxs) {
  struct MassMomsCalc {
    void operator()(const viewd_constgbx d_gbxs, Buffer<uint32_t>::mirrorviewd_buffer d_mom0,
                    Buffer<float>::mirrorviewd_buffer d_mom1,
                    Buffer<float>::mirrorviewd_buffer d_mom2) const {
      calculate_massmoments(d_gbxs, d_mom0, d_mom1, d_mom2);
    }
  } calc;

  const auto names = std::array<std::string_view, 3>({"massmom0", "massmom1", "massmom2"});

  return ConstTstepObserver(interval, DoMassMomsObs(dataset, maxchunk, ngbxs, names, calc));
}

#endif  // LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_
