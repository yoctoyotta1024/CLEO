# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# essential CLEO libaries #
add_subdirectory(gridboxes)
add_subdirectory(initialise)
add_subdirectory(runcleo)
add_subdirectory(superdrops)
add_subdirectory(zarr)
add_subdirectory(observers)
add_subdirectory(zarr2)

# libraries particular to certain CLEO configurations #
add_subdirectory(coupldyn_cvode EXCLUDE_FROM_ALL)
add_subdirectory(coupldyn_fromfile EXCLUDE_FROM_ALL)
add_subdirectory(coupldyn_yac EXCLUDE_FROM_ALL)
add_subdirectory(cartesiandomain EXCLUDE_FROM_ALL)
