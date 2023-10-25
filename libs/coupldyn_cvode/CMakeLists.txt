# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  #cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler 
endif()

# set library name
set(LIBNAME "coupldyn_cvode")