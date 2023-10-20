/*
 * ----- CLEO -----
 * File: condensation.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 20th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct for modelling condensation
 * microphysical process in SDM
 */

#ifndef CONDENSATION_HPP
#define CONDENSATION_HPP

#include <iostream>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "./microphysicalprocess.hpp"

struct DoCondensation
{
private:
  KOKKOS_INLINE_FUNCTION
  void do_condensation(const unsigned int subt) const
  {
    std::cout << "cond microphys @ t = " << subt << "\n";
  }

public:

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int subt) const
  /* this operator is used as an "adaptor" for using
  condensation as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    do_condensation(subt);
  }

};

MicrophysicalProcess auto
Condensation(const unsigned int interval)
/* constructs SdmProcess for condensation with constant timestep 'interval'
given a function to convert the interval to a (dimensionless) time
and the arguments required to construct the condensation method */
{
  return ConstTstepMicrophysics(interval, DoCondensation{});
}

#endif // CONDENSATION_HPP