# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# fetch and make Kokkos library available
set(FLAGS_B4 "${CMAKE_CXX_FLAGS}")
include(FetchContent)
FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.2.00.tar.gz
  GIT_TAG 71a9bcae52543bd065522bf3e41b5bfa467d8015
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC") # Kokkos needs -fPIC flag
FetchContent_MakeAvailable(Kokkos)
set(CMAKE_CXX_FLAGS "${FLAGS_B4}")

message(STATUS "Kokkos installation in: ${CMAKE_CURRENT_BINARY_DIR}")

# print information about compilers for C++ (and possibly CUDA)
if(${Kokkos_ENABLE_CUDA})
  message(STATUS "Using Kokkos nvcc wrapper (see: \
https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Compiling.html?highlight=wrapper#building-for-cuda)")
  message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
  message(STATUS "CC compiler: ${CMAKE_C_COMPILER}")
  message(STATUS "wrapper default (C++) compiler: ${NVCC_WRAPPER_DEFAULT_COMPILER}")
  message(STATUS "wrapper CUDA compiler: ${CUDA_ROOT}/bin/nvcc")
else()
  message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
  message(STATUS "CC compiler: ${CMAKE_C_COMPILER}")
endif()
