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
#include <concepts>

#include "./observers.hpp"
#include "./write_gridboxes.hpp"

template <typename Store>
struct MassMomArrays {
  XarrayZarrArray<Store, uint32_t> mom0_xzarr;  ///< 0th mass moment array
  XarrayZarrArray<Store, float> mom1_xzarr;     ///< 1st mass moment array
  XarrayZarrArray<Store, float> mom2_xzarr;     ///< 2nd mass moment array

  Buffer<uint32_t>::viewh_buffer h_mom0;  // view on host for 0th mass moment from every gridbox
  Buffer<uint32_t>::mirrorviewd_buffer d_mom0;  // mirror view of h_mom0 on device
  Buffer<float>::viewh_buffer h_mom1;        // view on host for 1st mass moment from every gridbox
  Buffer<float>::mirrorviewd_buffer d_mom1;  // mirror view of h_mom1 on device
  Buffer<float>::viewh_buffer h_mom2;        // view on host for 2nd mass moment from every gridbox
  Buffer<float>::mirrorviewd_buffer d_mom2;  // mirror view of h_mom2 on device

  /* create array for 0th mass moment for 4 byte unsigned integers (uint32_t type) */
  inline XarrayZarrArray<Store, uint32_t> create_array(Dataset<Store> &dataset,
                                                       const size_t maxchunk,
                                                       const size_t ngbxs) const {
    const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
    return dataset.template create_array<uint32_t>("massmom0", "", "<u4", 1, chunkshape,
                                                   {"time", "gbxindex"});
  }

  /* create array for >0th mass moment for 4 byte floating point numbers. Note conversion of
  scale factor from double (8 bytes) to single precision (4 bytes float) */
  inline XarrayZarrArray<Store, float> create_array(
      Dataset<Store> &dataset, const std::string_view name, const std::string_view units,
      const double scale_factor, const size_t maxchunk, const size_t ngbxs) const {
    const auto chunkshape = good2Dchunkshape(maxchunk, ngbxs);
    const auto scale_factor_ = static_cast<float>(scale_factor);
    return dataset.template create_array<float>(name, units, "<f4", scale_factor_, chunkshape,
                                                {"time", "gbxindex"});
  }

  MassMomArrays(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs)
      : mom0_xzarr(create_array(dataset, maxchunk, ngbxs)),
        mom1_xzarr(create_array(dataset, "massmom1", "g", dlc::MASS0grams, maxchunk, ngbxs)),
        mom2_xzarr(create_array(dataset, "massmom2", "g^2", dlc::MASS0grams * dlc::MASS0grams,
                                maxchunk, ngbxs)),
        h_mom0("h_mom0", ngbxs),
        d_mom0(Kokkos::create_mirror_view(ExecSpace(), h_mom0)),
        h_mom1("h_mom1", ngbxs),
        d_mom1(Kokkos::create_mirror_view(ExecSpace(), h_mom1)),
        h_mom2("h_mom2", ngbxs),
        d_mom2(Kokkos::create_mirror_view(ExecSpace(), h_mom2)) {}

  void write_arrayshape(Dataset<Store> &dataset) {
    dataset.write_arrayshape(mom0_xzarr);
    dataset.write_arrayshape(mom1_xzarr);
    dataset.write_arrayshape(mom2_xzarr);
  }

  void write_massmoments_to_arrays(Dataset<Store> &dataset) {
    Kokkos::deep_copy(h_mom0, d_mom0);
    dataset.write_to_array(mom0_xzarr, h_mom0);

    Kokkos::deep_copy(h_mom1, d_mom1);
    dataset.write_to_array(mom1_xzarr, h_mom1);

    Kokkos::deep_copy(h_mom2, d_mom2);
    dataset.write_to_array(mom2_xzarr, h_mom2);
  }
};

/* template class for observing 0th, 1st and 2nd mass moment (i.e. 0th, 3rd and 6th radius moment)
of droplet distribution in each gridbox to arrays in a dataset in a store */
template <typename Store>
class DoMassMomsObs {
 private:
  Dataset<Store> &dataset;                           ///< dataset to write moments to
  std::shared_ptr<MassMomArrays<Store>> xzarrs_ptr;  ///< pointer to mass moment arrays in dataset

 public:
  DoMassMomsObs(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs)
      : dataset(dataset),
        xzarrs_ptr(std::make_shared<MassMomArrays<Store>>(dataset, maxchunk, ngbxs)) {}

  ~DoMassMomsObs() { xzarrs_ptr->write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes time observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    calculate_massmoments(d_gbxs);
    xzarrs_ptr->write_massmoments_to_arrays(dataset);
  }
};

/* constructs observer which writes mass moments of droplet distribution in each gridbox
with a constant timestep 'interval' using an instance of the ConstTstepObserver class */
template <typename Store>
inline Observer auto MassMomentsObserver(const unsigned int interval, Dataset<Store> &dataset,
                                         const int maxchunk, const size_t ngbxs) {
  return ConstTstepObserver(interval, DoMassMomsObs(dataset, maxchunk, ngbxs));
}

#endif  // LIBS_OBSERVERS2_MASSMOMENTS_OBSERVER_HPP_
