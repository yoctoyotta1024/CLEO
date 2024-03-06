cmake_minimum_required(VERSION 3.21.1)

set(FLAGS_B4 "${CMAKE_CXX_FLAGS}")

# fetch and make Kokkos library available
include(FetchContent)
FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.2.00.tar.gz
  GIT_TAG 71a9bcae52543bd065522bf3e41b5bfa467d8015
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
FetchContent_MakeAvailable(Kokkos)
message(STATUS "Kokkos installation in: ${CMAKE_BINARY_DIR}/_deps/Kokkos")

set(CMAKE_CXX_FLAGS "${FLAGS_B4}")
