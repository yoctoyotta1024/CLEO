/*
 * ----- CLEO -----
 * File: nulldynamics.hpp
 * Project: coupldyn_null
 * Created Date: Friday 17th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Friday 17th November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * struct obeying coupleddynamics concept
 * for dynamics solver in CLEO where there
 * is no dynamics nor coupling to SDM
 */

#ifndef NULLDYNAMICS_HPP
#define NULLDYNAMICS_HPP

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"

struct NullDynamics
{
private:
  const unsigned int interval;

public:
  NullDynamics(const unsigned int couplstep)
      : interval(couplstep) {}

  KOKKOS_INLINE_FUNCTION
  void prepare_to_timestep() const {}

  KOKKOS_INLINE_FUNCTION
  auto get_couplstep() const
  {
    return interval;
  }

  KOKKOS_INLINE_FUNCTION
  bool on_step(const unsigned int t_mdl) const { return false; }

  KOKKOS_INLINE_FUNCTION
  void run_step(const unsigned int t_mdl,
                const unsigned int t_next) const {}
};

#endif // NULLDYNAMICS_HPP