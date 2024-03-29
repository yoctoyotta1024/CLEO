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
  XarrayZarrArray<Store, size_t> m0_xzarr;  ///< 0th mass moment array
  XarrayZarrArray<Store, float> m1_xzarr;   ///< 1st mass moment array
  XarrayZarrArray<Store, float> m2_xzarr;   ///< 2nd mass moment array

  template <typename T>
  XarrayZarrArray<Store, T> create_array(Dataset<Store> &dataset, const std::string_view name,
                                         const std::string_view units, const std::string_view dtype,
                                         const double scale_factor,
                                         const std::vector<size_t> &chunkshape) {
    return dataset.template create_array<T>(name, units, dtype, scale_factor, chunkshape,
                                            {"time", "gbxindex"});
  }

  MassMomArrays(Dataset<Store> &dataset, std::vector<size_t> &chunkshape)
      : m0_xzarr(create_array<uint32_t>(dataset, "massmom0", "", "<u4", 1, chunkshape)),
        m1_xzarr(create_array<float>(dataset, "massmom1", "g", "<f4", dlc::MASS0grams, chunkshape)),
        m2_xzarr(create_array<float>(dataset, "massmom2", "g^2", "<f4",
                                     dlc::MASS0grams * dlc::MASS0grams, chunkshape)) {}

  void write_arrayshape(Dataset<Store> &dataset) const {
    dataset.write_arrayshape(m0_xzarr);
    dataset.write_arrayshape(m1_xzarr);
    dataset.write_arrayshape(m2_xzarr);
  }
}

/* template class for observing 0th, 1st and 2nd mass moment (i.e. 0th, 3rd and 6th radius moment)
of droplet distribution in each gridbox to arrays in a dataset in a store */
template <typename Store>
class DoMassMomsObs {
 private:
  Dataset<Store> &dataset;                           ///< dataset to write moments to
  std::shared_ptr<MassMomArrays<Store>> xzarrs_ptr;  ///< pointer to mass moment array in dataset

 public:
  DoMassMomsObs(Dataset<Store> &dataset, const size_t maxchunk, const size_t ngbxs)
      : dataset(dataset),
        xzarr_ptr(
            std::make_shared<MassMomArrays<Store>>(dataset, good2Dchunkshape(maxchunk, ngbxs))) {}

  ~DoMassMomsObs() { xzarrs_ptr.write_arrayshape(dataset); }

  void before_timestepping(const viewd_constgbx d_gbxs) const {
    std::cout << "observer includes time observer\n";
  }

  void after_timestepping() const {}

  void at_start_step(const unsigned int t_mdl, const viewd_constgbx d_gbxs,
                     const viewd_constsupers totsupers) const {
    // TODO(CB)
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
