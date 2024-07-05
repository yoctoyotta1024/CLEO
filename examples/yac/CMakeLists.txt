# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# subdirectories for specific yac examples of CLEO
add_subdirectory(fromfile/src EXCLUDE_FROM_ALL)
add_subdirectory(yac_3d/src EXCLUDE_FROM_ALL)
