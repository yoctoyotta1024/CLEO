# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  #cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler 
endif()

# subdirectories for specific examples of CLEO
add_subdirectory(adiabaticparcel/src EXCLUDE_FROM_ALL)
add_subdirectory(boxmodelcollisions/golovin/src EXCLUDE_FROM_ALL)
add_subdirectory(boxmodelcollisions/long/src EXCLUDE_FROM_ALL)
add_subdirectory(boxmodelcollisions/lowlist/src EXCLUDE_FROM_ALL)
add_subdirectory(constthermo2d/src EXCLUDE_FROM_ALL)
add_subdirectory(divfreemotion/src EXCLUDE_FROM_ALL)
add_subdirectory(rainshaft1d/src EXCLUDE_FROM_ALL)
add_subdirectory(speedtest/src EXCLUDE_FROM_ALL)