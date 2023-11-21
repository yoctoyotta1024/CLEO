/*
 * ----- CLEO -----
 * File: cvodedynamics.cpp
 * Project: coupldyn_cvode
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Tuesday 21st November 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality for coupleddynamics concept for
 * dynamics solver in CLEO where coupling is
 * two-way to cvode adiabatic parcel ODE solver
 */

#include "./cvodedynamics.hpp"

void CvodeDynamics::prepare_to_timestep() const
/* checks initial y has been set and then
prints statement about cvode ODEs configuration */
{
  if (y == NULL)
  {
    throw std::invalid_argument("Cvode y not initialised");
  }

  print_initODEstatement();
}

void CvodeDynamics::print_initODEstatement() const
/* print initial ODE setup to the terminal screen */
{
  const double dimless_next_t(step2dimlesstime(interval));

  std::cout << "-------- CVODE ODE configuration ------------\n"
            << "No. Variables (NVARS) = " << NVARS << '\n'
            << "No. Equations (neq)   = " << neq << '\n'
            << "integer tstep         = " << interval << '\n'
            << "dimensionless tstep   = " << dimless_next_t << '\n'
            << "y0      = " << NV_Ith_S(y, 0) << '\n'
            << "y1      = " << NV_Ith_S(y, 1) << '\n'
            << "y2      = " << NV_Ith_S(y, 2) << '\n'
            << "y3      = " << NV_Ith_S(y, 3) << '\n'
            << "RTOL    = " << RTOL << '\n'
            << "ATOLS   = " << NV_Ith_S(ATOLS, 0) << '\n'
            << "---------------------------------------------\n\n";
};

int CvodeDynamics::run_dynamics(const unsigned int t_next)
/* make y before timestep new previousstate and then
integrate ODES for y from (dimensionless) time, t to
next_t = step2dimlesstime(t_next) */
{
  for (size_t i = 0; i < neq; ++i)
  {
    previousstates.at(i) = NV_Ith_S(y, i); // state 
  }

  const double dimless_next_t(step2dimlesstime(t_next));
  retval = CVode(cvode_mem, dimless_next_t, y, &t, CV_NORMAL);
  if (check_retval(&retval, "CVode", 1))
    return 1;

  return 0;
}

int CvodeDynamics::reinitialise(const double next_t,
                                const std::vector<double> &delta_y)
/* Reinitialize the solver after discontinuous change in
  temp, qv and qc due to condensation */
{
  re_y = NULL;
  re_y = N_VNew_Serial(neq, sunctx);

  for (size_t i = 0; i < neq; ++i)
  {
    NV_Ith_S(re_y, i) = NV_Ith_S(y, i) + delta_y[i];
  }

  retval = CVodeReInit(cvode_mem, next_t, re_y);
  if (check_retval((void *)&retval, "CVodeReInit", 1))
    return (1);

  return retval;
}

CvodeDynamics::CvodeDynamics(const Config &config,
                             const unsigned int couplstep,
                             const std::function<double(unsigned int)> step2dimlesstime)
    /* construct instance of CVODE ODE solver with initial conditions */
    : interval(couplstep),
      step2dimlesstime(step2dimlesstime),
      A(NULL),
      LS(NULL),
      cvode_mem(NULL),
      retval(0),
      neq(NVARS*config.ngbxs),
      t(0.0),
      y(NULL),
      re_y(NULL),
      ATOLS(NULL)
{
  data = (UserData)malloc(sizeof *data);
  previousstates = initial_conditions(config);

  const double wmax = (M_PI / 2) * (config.W_AVG / dlc::W0);  // dimensionless w velocity passed to thermo ODEs eg. dp_dt(t,y,ydot,w,...)
  const double tauhalf = (config.T_HALF / dlc::TIME0) / M_PI; // dimensionless timescale for w sinusoid
  init_userdata(neq, wmax, tauhalf);
  setup_ODE_solver(config.cvode_rtol, config.cvode_atol);
}

CvodeDynamics::~CvodeDynamics()
{
  /* print final statistics to the terminal screen */
  std::cout << "\nLast Iteration Statistics:\n";
  retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  /* free memory */
  N_VDestroy(y);            /* Free y vector */
  N_VDestroy(ATOLS);        /* Free abstol vector */
  free(data);               /* free user_data pointer struc */
  CVodeFree(&cvode_mem);    /* Free CVODE memory */
  SUNLinSolFree(LS);        /* Free the linear solver memory */
  SUNMatDestroy(A);         /* Free the matrix memory */
  SUNContext_Free(&sunctx); /* Free the SUNDIALS context */
}

std::vector<double>
CvodeDynamics::initial_conditions(const Config &config) const
/* return vector of dimensionless initial conditions
for thermodynamic variables (p, temp, qv, qc) to
initialise cvode thermodynamics solver */
{
  const double press_i(config.P_INIT / dlc::P0);
  const double temp_i(config.TEMP_INIT / dlc::TEMP0);
  const double qcond_i(0.0);

  const double psat(cvode_saturationpressure(temp_i));
  const double vapp(psat * config.relh_init / 100.0); // initial vapour pressure
  const double qvap_i(cvode_massmixingratio(vapp, press_i));

  std::vector<double> y_init(neq);
  for (size_t k = 0; k < neq; k += NVARS)
  {
    y_init.at(k) = press_i;
    y_init.at(k + 1) = temp_i;
    y_init.at(k + 2) = qvap_i;
    y_init.at(k + 3) = qcond_i;
  }

  return y_init;
}

void CvodeDynamics::init_userdata(const size_t neq,
                                  const double wmax,
                                  const double tauhalf)
/* set values in UserData structure for odes_func */
{
  data->neq = neq;
  data->wmax = wmax;
  data->tauhalf = tauhalf;
};

int CvodeDynamics::setup_ODE_solver(const double i_rtol,
                                    const double i_atol)
/* function does all the setup steps in order
to use CVODE sundials ODE solver */
{
  /* 0. Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1))
  {
    retval = 1;
  }

  /*  1. Initialize parallel or multi-threaded environment */
  // ------------------- (optional) --------------------- //

  /* 2. Set the scalar relative and vector absolute tolerances */
  RTOL = i_rtol;
  ATOLS = N_VNew_Serial(neq, sunctx);
  if (check_retval((void *)ATOLS, "N_VNew_Serial", 0))
    return (1);

  for (size_t i = 0; i < neq; ++i)
  {
    NV_Ith_S(ATOLS, i) = i_atol;
  }

  /* 3. initialise y vector with initial conditions */
  y = N_VNew_Serial(neq, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0))
    return (1);
  for (size_t i = 0; i < neq; ++i)
  {
    NV_Ith_S(y, i) = previousstates.at(i);
  }

  /* 4. Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula (CV_BDF) */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0))
    return (1);

  /* 5. Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0=0.0,
   * and the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, odes_func, 0.0, y);
  if (check_retval(&retval, "CVodeInit", 1))
    return (1);

  /* 6. Set linear solver optional inputs.
   * Provide user data which can be accessed in user provided routines */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval((void *)&retval, "CVodeSetUserData", 1))
    return (1);

  /* 7. Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, RTOL, ATOLS);
  if (check_retval(&retval, "CVodeSVtolerances", 1))
    return (1);

  /* 8. Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(neq, neq, sunctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0))
    return (1);

  /* 9. Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0))
    return (1);

  /* 10. Attach the matrix and linear solver to CVODE */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1))
    return (1);

  return 0;
};

int CvodeDynamics::check_retval(void *returnvalue,
                                const char *funcname,
                                int opt)
/* Check function return value...
  opt == 0 means SUNDIALS function allocates memory so check if
           returned NULL pointer
  opt == 1 means SUNDIALS function returns an integer value so check if
           retval < 0
  opt == 2 means function allocates memory so check if returned
           NULL pointer
*/
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    std::cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n"
              << funcname << " \n";
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int *)returnvalue;
    if (*retval < 0)
    {
      std::cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed with retval = %d\n\n"
                << funcname << " " << *retval << '\n';
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    std::cout << stderr << "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n"
              << funcname << '\n';
    return (1);
  }

  return (0);
}

