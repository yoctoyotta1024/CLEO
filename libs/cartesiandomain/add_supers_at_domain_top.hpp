/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: add_supers_at_domain_top.hpp
 * Project: cartesiandomain
 * Created Date: Tuesday 16th April 2024
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 20th April 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * A definition of the Domain Boundary Conditions to use for Cartesian GridBox Maps, Motion of
 * Super-Droplets and MoveSupersInDomain
 */

#ifndef LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
#define LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_

#include <Kokkos_Core.hpp>
#include <array>
#include <memory>
#include <stdexcept>

#include "../cleoconstants.hpp"
#include "../kokkosaliases.hpp"
#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/domainboundaries.hpp"
#include "gridboxes/sortsupers.hpp"
#include "initialise/optional_config_params.hpp"
#include "superdrops/superdrop.hpp"

struct AddSupersAtDomainTop {
 private:
  size_t newnsupers; /**< number of superdroplets to add to gridboxes above coord3lim */
  double coord3lim;  /**< gridboxes with upper bound > coord3lim get new super-droplets */
  std::shared_ptr<Superdrop::IDType::Gen>
      sdIdGen; /**< Pointer Superdrop::IDType object for super-droplet ID generation. */

  /* set super-droplet sdgbxindex to out of bounds value */
  void remove_supers_from_gridbox(const Gridbox &gbx) const;

  void add_supers_for_gridbox(const CartesianMaps &gbxmaps, const Gridbox &gbx,
                              const viewd_supers totsupers) const;

  Superdrop create_superdrop(const Gridbox &gbx) const;

  std::array<double, 3> create_superdrop_coords() const;

  SuperdropAttrs create_superdrop_attrs() const;

  /* (re)sorting supers based on their gbxindexes and then updating the span for each
  gridbox accordingly.
  Kokkos::parallel_for([...]) (on host) is equivalent to:
  for (size_t ii(0); ii < ngbxs; ++ii){[...]}
  when in serial */
  void move_supers_between_gridboxes(const viewd_gbx d_gbxs, const viewd_supers totsupers) const;

 public:
  /* New super-droplets are added to domain with coord3 >= COORD3LIM [m]. Note generation of
   * nextsdId assumes it is the only method creating super-droplets - otherwise created sdId may not
   * be unique*/
  explicit AddSupersAtDomainTop(const OptionalConfigParams::AddSupersAtDomainTopParams &config)
      : newnsupers(config.newnsupers),
        coord3lim(config.COORD3LIM / dlc::COORD0),
        sdIdGen(std::make_shared<Superdrop::IDType::Gen>(config.initnsupers)) {}

  /*_Note:_ totsupers is view of all superdrops (both in and out of bounds of domain).*/
  void operator()(const CartesianMaps &gbxmaps, viewd_gbx d_gbxs,
                  const viewd_supers totsupers) const;
};

#endif  // LIBS_CARTESIANDOMAIN_ADD_SUPERS_AT_DOMAIN_TOP_HPP_
