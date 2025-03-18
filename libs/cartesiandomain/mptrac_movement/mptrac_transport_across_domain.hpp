/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: mptrac_transport_across_domain.cpp
 * Project: mptrac_movement
 * Created Date: Monday 24th Febuary 2025
 * Author: Clara Bayley (CB)
 * Additional Contributors: Jan Clemens (JC)
 * -----
 * Last Modified: Monday 24th Febuary 2025
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Movement of a superdroplet throughout a cartesian domain using the MPTRAC library
 */

#ifndef LIBS_CARTESIANDOMAIN_MPTRAC_MOVEMENT_MPTRAC_TRANSPORT_ACROSS_DOMAIN_HPP_
#define LIBS_CARTESIANDOMAIN_MPTRAC_MOVEMENT_MPTRAC_TRANSPORT_ACROSS_DOMAIN_HPP_

#include <Kokkos_Core.hpp>
#include <array>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <memory>
#include <stdexcept>

#include "../../cleoconstants.hpp"
#include "../../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "gridboxes/gridbox.hpp"
#include "gridboxes/gridboxmaps.hpp"
#include "gridboxes/supersindomain.hpp"
#include "mpi.h"
#include "superdrops/superdrop.hpp"

struct MPTRACTransportAcrossDomain {
  /* define MPI_Datatype for MPTRAC superdroplet particle */
  void register_MPI_type_particle(MPI_Datatype *MPI_Particle, const int nquantities) const;

  /*
  function to move super-droplets between MPI processes using MPTRAC library,
  e.g. for superdroplets which move to/from gridboxes on different nodes.
  */
  viewd_supers sendrecv_supers(const CartesianMaps &gbxmaps, const viewd_gbx d_gbxs,
                              viewd_supers totsupers) const;

  /* (re)sorting supers based on their gbxindexes as step to 'move' superdroplets across the domain.
  May also include MPI communication with moves superdroplets away from/into a node's domain
  */
  SupersInDomain operator()(const CartesianMaps &gbxmaps, const viewd_gbx d_gbxs,
                            SupersInDomain &allsupers) const;
};

#endif  // LIBS_CARTESIANDOMAIN_MPTRAC_MOVEMENT_MPTRAC_TRANSPORT_ACROSS_DOMAIN_HPP_
