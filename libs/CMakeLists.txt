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
add_subdirectory(coupled_dynamics)

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
