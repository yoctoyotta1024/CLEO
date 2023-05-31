// Author: Clara Bayley
// File: cvodethermosolver.cpp
/* This file contains the CVODE ode solver
functionality for the evolution of the thermodynamics
(p, temp, qv and qc) over time */

#include "cvodethermosolver.hpp"

CvodeThermoSolver::CvodeThermoSolver()
    /* construct instance of CVODE ODE solver with
    initialised NULL vectors, matrix and solver */
    : A(NULL), LS(NULL), cvode_mem(NULL), retval(0),
    neq(NVARS), t(0.0), y(NULL), re_y(NULL), ATOLS(NULL)
{
  data = (UserData)malloc(sizeof *data);
}

CvodeThermoSolver::CvodeThermoSolver(const Config &config,
                                     const std::vector<double> &yinit_vector)
    /* construct instance of CVODE ODE
    solver with initial conditions */
    : A(NULL), LS(NULL), cvode_mem(NULL), retval(0), 
    neq(yinit_vector.size()), t(0.0), y(NULL), re_y(NULL), ATOLS(NULL)
{
  data = (UserData)malloc(sizeof *data);

  const double wmax = (M_PI / 2) * (config.W_AVG / dlc::W0);  // dimensionless w velocity passed to thermo ODEs eg. dp_dt(t,y,ydot,w,...)
  const double tauhalf = (config.T_HALF / dlc::TIME0) / M_PI; // dimensionless timescale for w sinusoid
  init_userdata(neq, config.doThermo, wmax, tauhalf);

  const std::array<double, NVARS> var_atols = {config.cvode_atol_p,
                                            config.cvode_atol_temp,
                                            config.cvode_atol_qv,
                                            config.cvode_atol_qc};
  std::vector<double> atols(neq);
  for (size_t k = 0; k < neq; k += NVARS)
  {
    atols.at(k) =  std::get<0>(var_atols);
    atols.at(k + 1) = std::get<1>(var_atols);
    atols.at(k + 2) = std::get<2>(var_atols);
    atols.at(k + 3) = std::get<3>(var_atols);
  }

  setup_ODE_solver(config.cvode_rtol, atols, yinit_vector);
}

CvodeThermoSolver::~CvodeThermoSolver()
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

void CvodeThermoSolver::init_userdata(const size_t neq,
                                      const bool doThermo,
                                      const double wmax,
                                      const double tauhalf)
/* set values in UserData structure for odes_func */
{
  data->neq = neq;
  data->doThermo = doThermo;
  data->wmax = wmax;
  data->tauhalf = tauhalf;
};

int CvodeThermoSolver::setup_ODE_solver(const double i_rtol,
                                        const std::vector<double> &i_atols,
                                        const std::vector<double> &i_yinit)
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
    NV_Ith_S(ATOLS, i) = i_atols.at(i);
  }

  /* 3. initialise y vector with initial conditions */
  y = N_VNew_Serial(neq, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0))
    return (1);
  for (size_t i = 0; i < neq; ++i)
  {
    NV_Ith_S(y, i) = i_yinit.at(i);
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

int CvodeThermoSolver::
    run_cvodestep(const int t_mdl, const int couplstep,
                  const double next_t)
/* Advance ODE solution in time to (dimless) next_t if on coupl step */
{
  if (t_mdl % couplstep == 0)
  {
    retval = CVode(cvode_mem, next_t, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1))
      return 1;
  }

  return 0;
}

int CvodeThermoSolver::reinitialise(const double next_t,
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

int CvodeThermoSolver::check_retval(void *returnvalue, const char *funcname, int opt)
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

int CvodeThermoSolver::print_init_ODEdata(const double tstep, const double tend)
/* print initial ODE setup to the terminal screen */
{

  if (y == NULL)
  {
    retval = -1;
    if (check_retval(&retval, "print_init_ODEdata", 1))
      return (1);
    return retval;
  }

  std::cout << "-------- CVODE initial conditions ------------\n"
            << "No. Variables (NVARS) = " << NVARS << '\n'
            << "No. Equations (neq) = " << neq << '\n'
            << "y0      = " << NV_Ith_S(y, 0) << '\n'
            << "y1      = " << NV_Ith_S(y, 1) << '\n'
            << "y2      = " << NV_Ith_S(y, 2) << '\n'
            << "y3      = " << NV_Ith_S(y, 3) << '\n'
            << "---------------------------------------------\n"
            << "RTOL        = " << RTOL << '\n'
            << "ATOLS[0]    = " << NV_Ith_S(ATOLS, 0) << '\n'
            << "ATOLS[1]    = " << NV_Ith_S(ATOLS, 1) << '\n'
            << "ATOLS[2]    = " << NV_Ith_S(ATOLS, 2) << '\n'
            << "ATOLS[3]    = " << NV_Ith_S(ATOLS, 3) << '\n'
            << "---------------------------------------------\n"
            << "inital t    = 0.0 \n"
            << "cvode tstep      = " << tstep << '\n'
            << "last t = " << tend << '\n'
            << "nout = " << ceil(tend / tstep) << '\n'
            << "---------------------------------------------\n\n";

  retval = 0;
  if (check_retval(&retval, "print_init_ODEdata", 1))
    return (1);
  return retval;
};
