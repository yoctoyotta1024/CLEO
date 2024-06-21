# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# essential CLEO libaries #
add_subdirectory(scratch EXCLUDE_FROM_ALL)
