# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# CLEO's CoupledDynamics Options
set(valid_coupled_dynamics "null" "fromfile" "cvode" "yac" "numpy")
if("${CLEO_COUPLED_DYNAMICS}" IN_LIST valid_coupled_dynamics)
  if("${CLEO_COUPLED_DYNAMICS}" STREQUAL "null")
    message(STATUS "CLEO using no coupled dynamics (${CLEO_COUPLED_DYNAMICS})")
  elseif("${CLEO_COUPLED_DYNAMICS}" STREQUAL "fromfile")
    add_subdirectory(coupldyn_fromfile EXCLUDE_FROM_ALL)
  elseif("${CLEO_COUPLED_DYNAMICS}" STREQUAL "cvode")
    add_subdirectory(coupldyn_cvode EXCLUDE_FROM_ALL)
  elseif("${CLEO_COUPLED_DYNAMICS}" STREQUAL "yac")
    add_subdirectory(coupldyn_yac EXCLUDE_FROM_ALL)
  # "numpy" case contained within add_subdirectory(cleo_python_bindings [...]), see below
  endif()
elseif("${CLEO_COUPLED_DYNAMICS}" STREQUAL "all") # all
  add_subdirectory(coupldyn_cvode EXCLUDE_FROM_ALL)
  add_subdirectory(coupldyn_fromfile EXCLUDE_FROM_ALL)
  add_subdirectory(coupldyn_yac EXCLUDE_FROM_ALL)
  # "numpy" case also included if add_subdirectory(cleo_python_bindings [...]), see below
elseif("${CLEO_COUPLED_DYNAMICS}" STREQUAL "") # default
  add_subdirectory(coupldyn_cvode EXCLUDE_FROM_ALL)
  add_subdirectory(coupldyn_fromfile EXCLUDE_FROM_ALL)
else()
  set(errmsg "Invalid CLEO_COUPLED_DYNAMICS '${CLEO_COUPLED_DYNAMICS}', please choose from: ${valid_coupled_dynamics}")
  message(FATAL_ERROR "${errmsg}")
endif()
