/*
 * ----- CLEO -----
 * File: createsupers.cpp
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


#include "./createsupers.hpp"

viewd_supers CreateSupers::operator()() const
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

viewd_supers CreateSupers::initialise_supers() const
/* initialise "totnsupers" number of superdrops in
kokkos view on device using InitSupers instance
for their initial gbxindex, spatial coordinates
and attributes */
{
  const size_t totnsupers(is.get_totnsupers());
  const InitSupersData isd(is);
  auto sdIdGen = Superdrop::IDType::Gen{};
  
  viewd_supers supers("supers", totnsupers);
  for (size_t kk(0); kk < totnsupers; ++kk)
  {
    const unsigned int sd_gbxindex(isd.sd_gbxindex.at(kk));
    const double coord3(isd.coord3.at(kk));
    const double coord1(isd.coord1.at(kk));
    const double coord2(isd.coord2.at(kk));
    const SuperdropAttrs attrs(isd.attrs_at(kk));
    const auto sd_id = sdIdGen.next();

    supers(kk) = Superdrop(sd_gbxindex, coord3,
                           coord1, coord2,
                           attrs, sd_id);
  }

  ensure_initialisation_complete(supers); 
  
  return supers;
}

void CreateSupers::
    ensure_initialisation_complete(viewd_constsupers supers) const
/* create "totnsupers" number of superdrops in view on device
using initsupers instance for initial data */
{
  if (supers.extent(0) < is.get_size())
  {
    const std::string err("Fewer superdroplets were created than were"
                          " given by initialisation data ie. " +
                          std::to_string(supers.extent(0)) + " < " +
                          std::to_string(is.get_size()));
    throw std::invalid_argument(err);
  }
}

void CreateSupers::print_supers(viewd_constsupers supers) const
{
  for (size_t kk(0); kk < supers.extent(0); ++kk)
  { 
    std::cout << "sdid: " << supers(kk).id.value << "\n"; 
  }
}

viewd_supers CreateSupers::sort_supers(viewd_supers supers) const
{
  return supers;
}