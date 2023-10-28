/*
 * ----- CLEO -----
 * File: differentialfuncs.hpp
 * Project: coupldyn_cvode
 * Created Date: Saturday 28th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Saturday 28th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
/* Header file for ODEs which are solved by
 CVODE ode solver to model evolution of the
thermodynamics (p, temp, qv and qc) over time */

#ifndef DIFFERENTIALFUNCS_HPP
#define DIFFERENTIALFUNCS_HPP

#include <cmath>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */

#include "cleoconstants.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/* user data structure for passing
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct pUserData
{
  size_t neq;
  realtype wmax;
  realtype tauhalf;
} * UserData;

int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data);
/* Simple function f(t,y, ydot) called by ODE solver to
  integrate ODEs over time. */


#endif // DIFFERENTIALFUNCS_HPP


