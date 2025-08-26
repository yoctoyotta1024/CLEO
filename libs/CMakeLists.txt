# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# essential CLEO libaries #
add_subdirectory(configuration)
add_subdirectory(initialise)
add_subdirectory(superdrops)
add_subdirectory(gridboxes)
add_subdirectory(runcleo)
add_subdirectory(zarr)
add_subdirectory(observers)

# libraries particular to certain CLEO configurations #
# CLEO's CoupledDynamics
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

# CLEO's Domain
set(valid_domain "cartesian")
if("${CLEO_DOMAIN}" IN_LIST valid_domain)
  if("${CLEO_DOMAIN}" STREQUAL "cartesian")
    add_subdirectory(cartesiandomain EXCLUDE_FROM_ALL)
  endif()
elseif("${CLEO_DOMAIN}" STREQUAL "" OR "${CLEO_DOMAIN}" STREQUAL "all") # default and all
  add_subdirectory(cartesiandomain EXCLUDE_FROM_ALL)
else()
  set(errmsg "Invalid CLEO_DOMAIN '${CLEO_DOMAIN}', please choose from: ${valid_domain}")
  message(FATAL_ERROR "${errmsg}")
endif()

# optionally make CLEO's python bindings
if(CLEO_NO_PYBINDINGS)
  message(STATUS "CLEO excluding python bindings CLEO_NO_PYBINDINGS=${CLEO_NO_PYBINDINGS}")
else()
  add_subdirectory(cleo_python_bindings EXCLUDE_FROM_ALL)
endif()
