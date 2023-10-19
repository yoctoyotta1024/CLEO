# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  #cmake_minimum_required(VERSION 3.21.1) # if using Kokkos with NVC++ compiler 
endif()

# essential CLEO SDM libaries #
add_subdirectory(initialise)
add_subdirectory(observers)
add_subdirectory(runcleo)
add_subdirectory(gridboxes)
add_subdirectory(superdrops)
add_subdirectory(zarr)

# libraries particular to certain CLEO configurations #
add_subdirectory(coupldyn_cvode EXCLUDE_FROM_ALL)
add_subdirectory(coupldyn_fromfile EXCLUDE_FROM_ALL)
add_subdirectory(cartesiandomain EXCLUDE_FROM_ALL)