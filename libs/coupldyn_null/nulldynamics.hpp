/* Copyright (c) 2023 MPI-M, Clara Bayley
 *
 * ----- CLEO -----
 * File: nulldynamics.hpp
 * Project: coupldyn_null
 * Created Date: Friday 17th November 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * struct obeying coupleddynamics concept
 * for dynamics solver in CLEO where there
 * is no dynamics nor coupling to SDM
 */

#ifndef LIBS_COUPLDYN_NULL_NULLDYNAMICS_HPP_
#define LIBS_COUPLDYN_NULL_NULLDYNAMICS_HPP_

#include <Kokkos_Core.hpp>

#include "../kokkosaliases.hpp"

struct NullDynamics {
 private:
  const unsigned int interval;

 public:
  explicit NullDynamics(const unsigned int couplstep) : interval(couplstep) {}

  void prepare_to_timestep() const {}

  unsigned int get_couplstep() const { return interval; }

  bool on_step(const unsigned int t_mdl) const { return false; }

  void run_step(const unsigned int t_mdl, const unsigned int t_next) const {}
};

#endif  // LIBS_COUPLDYN_NULL_NULLDYNAMICS_HPP_
