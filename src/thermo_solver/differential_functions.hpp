// Author: Clara Bayley
// File: differential_functions.hpp
/* Header file for ODEs which are solved by
 CVODE ode solver to model evolution of the
thermodynamics (p, temp, qv and qc) over time */

#ifndef DIFFERENTIAL_FUNCTIONS_HPP
#define DIFFERENTIAL_FUNCTIONS_HPP

// #include <iostream>
#include <cmath>
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */

#include "../claras_SDconstants.hpp"

namespace dlc = dimless_constants;
namespace DC = dimmed_constants;

/* user data structure for passing
      args to f() function from ode solver */
/* Type : UserData contains preconditioner blocks,
     pivot arrays, and problem constants */
typedef struct pUserData
{
  size_t neq;
  bool doThermo;
  realtype wmax;
  realtype tauhalf;
} * UserData;

int odes_func(realtype t, N_Vector y, N_Vector ydot, void *user_data);
/* Simple function f(t,y, ydot) called by ODE solver to
  integrate ODEs over time. */

#endif // DIFFERENTIAL_FUNCTIONS_HPP