# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  # cmake_minimum_required(VERSION 3.18.0)
  cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler 
endif()

# set project name and source directory of .cpp files
project("adiabatic_parcel")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# set(compileops -Wall -pedantic -g -gdwarf-4) # less strict without -Werror
# set compiler options
if(NOT DEFINED compileops)
  set(compileops ${defaultcompileops})
endif()

add_subdirectory(src)