// Author: Clara Bayley
// File: cvodethermosolver.hpp
/* Header file for (CVODE) ode solver
which models evolution of the thermodynamics
(p, temp, qv and qc) over time */

#ifndef CVODETHERMOSOLVER_HPP
#define CVODETHERMOSOLVER_HPP

#include <array>
#include <vector>
#include <algorithm>
#include <iostream>

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */

#include "initialisation/config.hpp"
#include "./differential_functions.hpp"

namespace dlc = dimless_constants;

class CvodeThermoSolver
{
private:
  /* SUNDIALS CVODE solver stuff */
  SUNContext sunctx;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval;

  /* ODE problem stuff */
  static constexpr int NVARS = 4;    // no. of distinct variables (= no. ODEs per grid box)
  const size_t neq;              // No. of equations/ODEs (= total no. variables across all Grid Boxes)
  realtype t;
  realtype RTOL;
  N_Vector y;
  N_Vector re_y;
  N_Vector ATOLS;
  UserData data;

  int check_retval(void *returnvalue, const char *funcname, int opt);
  /* Check function return value for memory or sundials CVODE error */

public:
  CvodeThermoSolver();

  CvodeThermoSolver(const Config &config, const std::vector<double> &yinit_vector);

  ~CvodeThermoSolver();
  /* print final statistics to the terminal
  screen and free CVODE memory */

  void init_userdata(const size_t neq, const bool doThermo,
                     const double wmax, const double tauhalf);
  /* set values in UserData structure for odes_func */

  int setup_ODE_solver(const double i_rtol,
                       const std::vector<double> &i_atols,
                       const std::vector<double> &i_yinit);
  /* function does all the setup steps in order
  to use CVODE sundials ODE solver */

  int run_cvodestep(const int t_mdl, const int couplstep,
                    const double next_t);
  /* Advance ODE solution in time to (dimless)
  next_t if on coupl step */

  int reinitialise(const double next_t,
                   const std::vector<double> &delta_y);
  /* Reinitialize the solver after discontinuous change in y*/

  int print_init_ODEdata(const double tstep, const double tend);
  /* print initial ODE setup to the terminal screen */

  inline double get_time() const { return t; }

  inline double get_pressure(const long unsigned int ii)
      const { return NV_Ith_S(y, NVARS * ii); }

  inline double get_temperature(const long unsigned int ii)
      const { return NV_Ith_S(y, NVARS * ii + 1); }

  inline double get_qvap(const long unsigned int ii)
      const { return NV_Ith_S(y, NVARS * ii + 2); }

  inline double get_qcond(const long unsigned int ii)
      const { return NV_Ith_S(y, NVARS * ii + 3); }
};

#endif // CVODETHERMOSOLVER_HPP