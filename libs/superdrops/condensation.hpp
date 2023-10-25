/*
 * ----- CLEO -----
 * File: condensation.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Wednesday 25th October 2023
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

#include <concepts>

#include <Kokkos_Core.hpp>

#include "./kokkosaliases_sd.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"
#include "./state.hpp"
#include "./urbg.hpp"

struct DoCondensation
/* function-like type that enacts
condensation / evaporation microphysical process */
{
private:
  KOKKOS_FUNCTION
  subviewd_supers do_condensation(const unsigned int subt,
                                  const subviewd_supers supers) const;

public:
  template <class DeviceType>
  KOKKOS_INLINE_FUNCTION
  subviewd_supers operator()(const unsigned int subt,
                             subviewd_supers supers,
                             State &state,
                             URBG<DeviceType> urbg) const
  /* this operator is used as an "adaptor" for using
  condensation as the MicrophysicsFunction type in a
  ConstTstepMicrophysics instance (*hint* which itself
  satsifies the MicrophysicalProcess concept) */
  {
    supers = do_condensation(subt, supers);

    return supers;
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