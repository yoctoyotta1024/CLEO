/*
 * ----- CLEO -----
 * File: massmomentobs.hpp
 * Project: observers
 * Created Date: Sunday 22nd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 22nd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * Observer to output nsupers per gridbox
 * to array in a zarr file system storage
 */

#ifndef MASSMOMENTOBS_HPP
#define MASSMOMENTOBS_HPP

#include <concepts>
#include <memory>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "./observers.hpp"
#include "superdrops/superdrop.hpp"
#include "gridboxes/gridbox.hpp"
#include "zarr/twodstorage.hpp"

inline Observer auto
MassMomentsObserver(const unsigned int interval,
                   FSStore &store,
                   const int maxchunk,
                   const size_t ngbxs);
/* constructs observer of the nth mass moment
'mom' in each gridbox with a constant
timestep 'interval' using an instance of the
DoMassMomentObs class */

class DoMassMomentsObs 
/* observe nth mass moment in each gridbox and
write it to an array 'zarr' store as determined
by the 2DStorage instance */
{
private:
  std::shared_ptr<MassMomentsStorage<double>> zarr;
  
public:
   DoMassMomentsObs(FSStore &store,
               const int maxchunk,
               const size_t ngbxs)
      : zarr(std::make_shared<store_type>(store, maxchunk,
                                         "<f8", ngbxs))
  {
    zarr->is_dim1(ngbxs, "gbxindex");
  }

  void before_timestepping(const viewh_constgbx h_gbxs) const
  {
    std::cout << "observer includes MassMomentsObserver\n";
  }

  void at_start_step(const unsigned int t_mdl,
                     const viewh_constgbx h_gbxs) const
  /* deep copy if necessary (if superdrops are on device not
  host memory), then writes mass moments to 2-D zarr storages */
  {
    const size_t ngbxs(h_gbxs.extent(0));
    for (size_t ii(0); ii < ngbxs; ++ii)
    {
      auto h_supers = h_gbxs(ii).hostcopy(); // deep copy if supers not on host memory
      massmoments_to_storage(h_supers);
    }
    ++(zarr->nobs);
  }

  void massmoments_to_storage(const mirrorh_constsupers h_supers)
  /* calculated 0th, 1st and 2nd moment of the (real) droplet mass
  distribution and then writes them to zarr storage. (I.e.
  0th, 3rd and 6th moment of the droplet radius distribution) */
  {
    double mom0(0.0); // 0th moment = number of (real) droplets
    double mom1(0.0); // 1st moment = mass of (real) droplets
    double mom2(0.0); // 2nd moment = mass^2 of (real) droplets
    for (size_t kk(0); kk < h_supers.extent(0); ++kk)
    {
      const(double) xi(h_supers(kk).get_xi()); // cast multiplicity from unsigned int to double
      const double mass(h_supers(kk).mass());
      mom0 += xi;
      mom1 += xi * mass;
      mom2 += xi * mass * mass;
    }
    
    zarr->mom0_to_storage(mom0);
    zarr->mom1_to_storage(mom1);
    zarr->mom2_to_storage(mom2);
  }

};


#endif MASSMOMENTOBS_HPP
