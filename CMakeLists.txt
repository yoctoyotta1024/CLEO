# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# e.g. could run cmake with these compiler options:
# cmake -S ./ -B ./build -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CC_COMPILER=gcc -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON

# set cmake version
cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) ### if using Kokkos with nvc++ compiler on Levante

project(project
  LANGUAGES CXX
  DESCRIPTION "CLEO by Clara Bayley and other developers"
  VERSION 1.0
)

# fetch and make Kokkos library available
set(CMAKE_CXX_STANDARD "20")
include(FetchContent)
FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.2.00.tar.gz
  GIT_TAG 71a9bcae52543bd065522bf3e41b5bfa467d8015
)
set(FLAGS_B4 "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
FetchContent_MakeAvailable(Kokkos)
set(CMAKE_CXX_FLAGS "${FLAGS_B4}")

if(${Kokkos_ENABLE_CUDA})
  message(STATUS "Using Kokkos nvcc wrapper (see: \
https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Compiling.html?highlight=wrapper#building-for-cuda)")
  set(CMAKE_CXX_COMPILER ${CLEO_NVCC_WRAPPER} CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CC_COMPILER ${CLEO_CC_COMPILER} CACHE STRING "CC compiler" FORCE)
  set(NVCC_WRAPPER_DEFAULT_COMPILER ${CLEO_CXX_COMPILER})
  set(CUDA_ROOT ${CLEO_CUDA_ROOT})
  message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
  message(STATUS "CC compiler: ${CMAKE_CC_COMPILER}")
  message(STATUS "wrapper default (CXX) compiler set to: ${NVCC_WRAPPER_DEFAULT_COMPILER}")
  message(STATUS "wrapper nvcc compiler: ${CUDA_ROOT}/bin/nvcc")
else()
  set(CMAKE_CXX_COMPILER ${CLEO_CXX_COMPILER} CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CC_COMPILER ${CLEO_CC_COMPILER} CACHE STRING "CC compiler" FORCE)
  message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
  message(STATUS "CC compiler: ${CMAKE_CC_COMPILER}")
endif()

# print default cmake compiler flags
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# subdirectories of CLEO libray and main program
add_subdirectory(libs)
add_subdirectory(src)
#add_subdirectory(roughpaper)

# subdirectories for specific examples of CLEO
add_subdirectory(examples EXCLUDE_FROM_ALL)
