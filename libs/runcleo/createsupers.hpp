/*
 * ----- CLEO -----
 * File: createsupers.hpp
 * Project: runcleo
 * Created Date: Tuesday 17th October 2023
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
 * file for functions used to create
 * superdroplets given their initial data
 */

#ifndef CREATESUPERS_HPP
#define CREATESUPERS_HPP

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"
#include "superdrops/superdrop.hpp"

template <typename FetchInitData>
class CreateSupers
{
private:
  const FetchInitData &fisd;

  struct InitData
  {
    std::vector<unsigned int> sd_gbxindex;
    std::vector<double> coord3;
    std::vector<double> coord1;
    std::vector<double> coord2;
    std::vector<double> radius;
    std::vector<double> msol;
    std::vector<unsigned long long> xi;
    std::shared_ptr<const SoluteProperties> solute;

    InitData(const FetchInitData &fisd)
        : sd_gbxindex(fisd.sd_gbxindex()),
          coord3(fisd.coord3()),
          coord1(fisd.coord1()),
          coord2(fisd.coord2()),
          radius(fisd.radius()),
          msol(fisd.msol()),
          xi(fisd.xi()),
          solute() {}

    SuperdropAttrs attrs_at(const unsigned int kk) const
    {
      return {radius.at(kk), msol.at(kk), xi.at(kk), solute};
    }
  };

  viewd_supers initialise_supers() const
  /* initialise "totnsupers" number of superdrops in
  kokkos view on device using InitSupers instance
  for their initial gbxindex, spatial coordinates
  and attributes */
  {
    const size_t totnsupers(fisd.get_totnsupers());
    const InitData idat(fisd);
    auto sdIdGen = Superdrop::IDType::Gen{};

    viewd_supers supers("supers", totnsupers);
    for (size_t kk(0); kk < totnsupers; ++kk)
    {
      const unsigned int sd_gbxindex(idat.sd_gbxindex.at(kk));
      const double coord3(idat.coord3.at(kk));
      const double coord1(idat.coord1.at(kk));
      const double coord2(idat.coord2.at(kk));
      const SuperdropAttrs attrs(idat.attrs_at(kk));
      const auto sd_id = sdIdGen.next();

      supers(kk) = Superdrop(sd_gbxindex, coord3,
                             coord1, coord2,
                             attrs, sd_id);
    }

    ensure_initialisation_complete(supers);

    return supers;
  }

  void ensure_initialisation_complete(viewd_constsupers supers) const
  /* create "totnsupers" number of superdrops in view on device
  using initsupers instance for initial data */
  {
    if (supers.extent(0) < fisd.get_size())
    {
      const std::string err("Fewer superdroplets were created than were"
                            " given by initialisation data ie. " +
                            std::to_string(supers.extent(0)) + " < " +
                            std::to_string(fisd.get_size()));
      throw std::invalid_argument(err);
    }
  }

  void print_supers(viewd_constsupers supers) const
  {
    for (size_t kk(0); kk < supers.extent(0); ++kk)
    {
      std::cout << "sdid: " << supers(kk).id.value << "\n";
    }
  }

  viewd_supers sort_supers(viewd_supers supers) const
  {
    return supers;
  }

public:
  CreateSupers(const FetchInitData &fisd) : fisd(fisd) {}

  viewd_supers operator()() const
  /* create view of "totnsupers" number of superdrops in a
  Kokkos view on the device memory ordered by their gridbox indexes
  and all with the same solute properties data and initial
  conditons from the referenced InitSupers struct */
  {
    viewd_supers supers(initialise_supers());

    supers = sort_supers(supers);

    print_supers(supers);

    return supers;
  }
};

#endif // CREATESUPERS_HPP