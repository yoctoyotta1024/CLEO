/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: cvodedynamics.cpp
 * Project: coupldyn_cvode
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * functionality for coupleddynamics concept for dynamics solver in CLEO where coupling is
 * two-way to cvode adiabatic parcel ODE solver
 */

#include "coupldyn_cvode/cvodedynamics.hpp"

/* checks initial y has been set and then
prints statement about cvode ODEs configuration */
void CvodeDynamics::prepare_to_timestep() const {
  if (y == NULL) {
    throw std::invalid_argument("Cvode y not initialised");
  }

  print_initODEstatement();
}

/* print initial ODE setup to the terminal screen */
void CvodeDynamics::print_initODEstatement() const {
  const auto dimless_next_t = double{step2dimlesstime(interval)};

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
}

/* make y before timestep new previousstate and then
integrate ODES for y from (dimensionless) time, t to
next_t = step2dimlesstime(t_next) */
int CvodeDynamics::run_dynamics(const unsigned int t_next) {
  for (size_t i = 0; i < neq; ++i) {
    previousstates.at(i) = NV_Ith_S(y, i);  // state
  }

  const auto dimless_next_t = double{step2dimlesstime(t_next)};
  retval = CVode(cvode_mem, dimless_next_t, y, &t, CV_NORMAL);
  if (check_retval(&retval, "CVode", 1)) return 1;

  return 0;
}

/* Reinitialize the solver after discontinuous change in
  temp, qv and qc due to condensation */
int CvodeDynamics::reinitialise(const double next_t, const std::vector<double> &delta_y) {
  re_y = NULL;
  re_y = N_VNew_Serial(neq, sunctx);

  for (size_t i = 0; i < neq; ++i) {
    NV_Ith_S(re_y, i) = NV_Ith_S(y, i) + delta_y[i];
  }

  retval = CVodeReInit(cvode_mem, next_t, re_y);
  if (check_retval(reinterpret_cast<void *>(&retval), "CVodeReInit", 1)) return (1);

  return retval;
}

/* construct instance of CVODE ODE solver with initial conditions */
CvodeDynamics::CvodeDynamics(const OptionalConfigParams::CvodeDynamicsParams &config,
                             const unsigned int couplstep,
                             const std::function<double(unsigned int)> step2dimlesstime)
    : interval(couplstep),
      step2dimlesstime(step2dimlesstime),
      A(NULL),
      LS(NULL),
      cvode_mem(NULL),
      retval(0),
      neq(NVARS * config.ngbxs),
      t(0.0),
      y(NULL),
      re_y(NULL),
      ATOLS(NULL) {
  data = (UserData)malloc(sizeof *data);
  previousstates = initial_conditions(config);

  const auto wmax = double{
      (M_PI / 2) *
      (config.W_avg /
       dlc::W0)};  // dimensionless w velocity passed to thermo ODEs eg. dp_dt(t,y,ydot,w,...)
  const auto tauhalf =
      double{(config.TAU_half / dlc::TIME0) / M_PI};  // dimensionless timescale for w sinusoid
  init_userdata(neq, wmax, tauhalf);
  setup_ODE_solver(config.rtol, config.atol);
}

/* print final statistics to the terminal screen */
CvodeDynamics::~CvodeDynamics() {
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

/* return vector of dimensionless initial conditions
for thermodynamic variables (p, temp, qv, qc) to
initialise cvode thermodynamics solver */
std::vector<double> CvodeDynamics::initial_conditions(
    const OptionalConfigParams::CvodeDynamicsParams &config) const {
  const auto press_i = double{config.P_init / dlc::P0};
  const auto temp_i = double{config.TEMP_init / dlc::TEMP0};
  const auto qcond_i = double{0.0};

  const auto psat = double{cvode_saturationpressure(temp_i)};
  const auto vapp = double{psat * config.relh_init / 100.0};  // initial vapour pressure
  const auto qvap_i = double{cvode_massmixingratio(vapp, press_i)};

  std::vector<double> y_init(neq);
  for (size_t k = 0; k < neq; k += NVARS) {
    y_init.at(k) = press_i;
    y_init.at(k + 1) = temp_i;
    y_init.at(k + 2) = qvap_i;
    y_init.at(k + 3) = qcond_i;
  }

  return y_init;
}

/* set values in UserData structure for odes_func */
void CvodeDynamics::init_userdata(const size_t neq, const double wmax, const double tauhalf) {
  data->neq = neq;
  data->wmax = wmax;
  data->tauhalf = tauhalf;
}

/* function does all the setup steps in order
to use CVODE sundials ODE solver */
int CvodeDynamics::setup_ODE_solver(const double i_rtol, const double i_atol) {
  /* 0. Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) {
    retval = 1;
  }

  /*  1. Initialize parallel or multi-threaded environment */
  // ------------------- (optional) --------------------- //

  /* 2. Set the scalar relative and vector absolute tolerances */
  RTOL = i_rtol;
  ATOLS = N_VNew_Serial(neq, sunctx);
  if (check_retval(reinterpret_cast<void *>(ATOLS), "N_VNew_Serial", 0)) return (1);

  for (size_t i = 0; i < neq; ++i) {
    NV_Ith_S(ATOLS, i) = i_atol;
  }

  /* 3. initialise y vector with initial conditions */
  y = N_VNew_Serial(neq, sunctx);
  if (check_retval(reinterpret_cast<void *>(y), "N_VNew_Serial", 0)) return (1);
  for (size_t i = 0; i < neq; ++i) {
    NV_Ith_S(y, i) = previousstates.at(i);
  }

  /* 4. Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula (CV_BDF) */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval(reinterpret_cast<void *>(cvode_mem), "CVodeCreate", 0)) return (1);

  /* 5. Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the initial time T0=0.0,
   * and the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, odes_func, 0.0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return (1);

  /* 6. Set linear solver optional inputs.
   * Provide user data which can be accessed in user provided routines */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(reinterpret_cast<void *>(&retval), "CVodeSetUserData", 1)) return (1);

  /* 7. Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, RTOL, ATOLS);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return (1);

  /* 8. Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(neq, neq, sunctx);
  if (check_retval(reinterpret_cast<void *>(A), "SUNDenseMatrix", 0)) return (1);

  /* 9. Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval(reinterpret_cast<void *>(LS), "SUNLinSol_Dense", 0)) return (1);

  /* 10. Attach the matrix and linear solver to CVODE */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return (1);

  return 0;
}

/* Check function return value...
  opt == 0 means SUNDIALS function allocates memory so check if
           returned NULL pointer
  opt == 1 means SUNDIALS function returns an integer value so check if
           retval < 0
  opt == 2 means function allocates memory so check if returned
           NULL pointer
*/
int CvodeDynamics::check_retval(void *returnvalue, const char *funcname, int opt) {
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    std::cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n"
              << funcname << " \n";
    return (1);

    // Check if retval < 0
  } else if (opt == 1) {
    retval = reinterpret_cast<int *>(returnvalue);
    if (*retval < 0) {
      std::cout << stderr << "\nCVODE_SUNDIALS_ERROR: %s() failed with retval = %d\n\n"
                << funcname << " " << *retval << '\n';
      return (1);
    }
    // Check if function returned NULL pointer - no memory allocated
  } else if (opt == 2 && returnvalue == NULL) {
    std::cout << stderr << "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n"
              << funcname << '\n';
    return (1);
  }

  return (0);
}
