/*
 * ----- CLEO -----
 * File: creategbxs.hpp
 * Project: runcleo
 * Created Date: Wednesday 18th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 18th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * file for structure to create a dualview of
 * gridboxes from using some initial conditions
 */


#ifndef CREATEGBXS_HPP
#define CREATEGBXS_HPP

#include <iostream>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include "../kokkosaliases.hpp"
#include "gridboxes/gridbox.hpp"

class CreateGbxs 
{
private:
  template <typename FetchInitData>
  inline dualview_gbx initialise_gbxs(const FetchInitData &fid) const
  /* initialise a view of superdrops (on device memory)
  using data from an InitData instance for their initial
  gbxindex, spatial coordinates and attributes */
  {
    const size_t ngbxs(fid.get_ngbxs());
    
    dualview_gbx gbxs("gbxs", ngbxs);
  }

  void ensure_initialisation_complete(dualview_gbx supers,
                                      const size_t size) const;


  inline void print_gbxs(dualview_gbx gbxs) const
  /* print gridboxes information */
  {
    std::cout << "some gbx info\n";
  }

public:
  template <typename FetchInitData>
  dualview_gbx operator()(const FetchInitData fid) const
  {

    std::cout << "\n--- create gridboxes ---"
              << "\ninitialising";
    dualview_gbx gbxs(initialise_gbxs(fid));

    std::cout << "\nmiddlestep";
    
    ensure_initialisation_complete(gbxs, fid.get_size());
    print_gbxs(gbxs);
    std::cout << "\n--- create gridboxes: success ---\n";

    return gbxs;
  }
};

#endif // CREATEGBXS_HPP