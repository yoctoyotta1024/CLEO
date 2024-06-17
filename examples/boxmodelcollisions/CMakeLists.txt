# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# subdirectories for specific examples of CLEO
add_subdirectory(golovin/src EXCLUDE_FROM_ALL)
add_subdirectory(long/src EXCLUDE_FROM_ALL)
add_subdirectory(lowlist/src EXCLUDE_FROM_ALL)
add_subdirectory(szakallurbich/src EXCLUDE_FROM_ALL)
add_subdirectory(testikstraub/src EXCLUDE_FROM_ALL)
