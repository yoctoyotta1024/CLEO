/*
 * ----- CLEO -----
 * File: createsupers.hpp
 * Project: runcleo
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 17th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * file for functions used to create
 * superdroplets given their initial data
 */


#ifndef CREATESUPERS_HPP
#define CREATESUPERS_HPP

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "initialise/initconds.hpp"
#include "superdrops/superdrop.hpp"

class CreateSupers
{
private:
  const InitSupers &is;

  struct InitSupersData
  {
    std::vector<unsigned int> sd_gbxindex;
    std::vector<double> coord3;
    std::vector<double> coord1;
    std::vector<double> coord2;
    std::vector<double> radius;
    std::vector<double> msol;
    std::vector<unsigned long long> xi;
    std::shared_ptr<const SoluteProperties> solute;

    InitSupersData(const InitSupers &is)
        : sd_gbxindex(is.sd_gbxindex()),
          coord3(is.coord3()),
          coord1(is.coord1()),
          coord2(is.coord2()),
          radius(is.radius()),
          msol(is.msol()),
          xi(is.xi()),
          solute() {}

    SuperdropAttrs attrs_at(const unsigned int kk) const
    {
      return {radius.at(kk), msol.at(kk), xi.at(kk), solute};
    }
  };

  viewd_supers initialise_supers() const;

  void ensure_initialisation_complete(viewd_constsupers supers) const;

  void print_supers(viewd_constsupers supers) const;

  viewd_supers sort_supers(viewd_supers supers) const;

public:
  CreateSupers(const InitSupers &is) : is(is) {}

  viewd_supers operator()() const;
};

#endif // CREATESUPERS_HPP