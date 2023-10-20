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
 * and evaporation of water
 * microphysical process in SDM
 */

#ifndef CONDENSATION_HPP
#define CONDENSATION_HPP

#include <iostream>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "./microphysicalprocess.hpp"

struct DoCondensation
/* function-like type that enacts
condensation / evaporation microphysical process */
{
private:
  KOKKOS_FUNCTION
  void do_condensation(const unsigned int subt) const;

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

inline MicrophysicalProcess auto
Condensation(const unsigned int interval)
/* constructs Microphysical Process for
condensation/evaporation of superdroplets with a
constant timestep 'interval' given the
"do_condensation" function-like type */
{
  return ConstTstepMicrophysics(interval, DoCondensation{});
}

#endif // CONDENSATION_HPP