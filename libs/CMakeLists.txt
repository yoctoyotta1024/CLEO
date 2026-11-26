# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# ${CLEO_SOURCE_DIR}/libs is a shared include path for every CLEO library. Libraries
# can link against ``cleo_libs_include`` instead of each declaring themselves.
add_library(cleo_libs_include INTERFACE)
target_include_directories(cleo_libs_include INTERFACE "${CLEO_SOURCE_DIR}/libs")

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
