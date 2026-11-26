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
add_subdirectory(grids)

# optionally make CLEO's python bindings
if(CLEO_NO_PYBINDINGS)
  message(STATUS "CLEO excluding python bindings CLEO_NO_PYBINDINGS=${CLEO_NO_PYBINDINGS}")
else()
  add_subdirectory(cleo_python_bindings EXCLUDE_FROM_ALL)
endif()
