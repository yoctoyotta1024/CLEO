# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# e.g. could run cmake with these compiler options:
# cmake -S ./ -B ./build -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CC_COMPILER=gcc -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON

# set cmake version
cmake_minimum_required(VERSION 3.18.0)
# cmake_minimum_required(VERSION 3.21.1) ### if using Kokkos with nvc++ compiler on Levante

project(project
  LANGUAGES CXX
  DESCRIPTION "CLEO by Clara Bayley and other developers"
  VERSION 1.0
)

# ensure C++ compiler uses certain settings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CXX_EXTENSIONS ON)

# install Kokkos for project
message(STATUS "Using Kokkos installation from: ${CMAKE_SOURCE_DIR}/extern/kokkos")
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/kokkos)

# print default compiler flags
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# add directories of CLEO libray and main program
add_subdirectory(libs)
add_subdirectory(src)

# add directories for specific examples of CLEO
add_subdirectory(examples EXCLUDE_FROM_ALL)

# add directory for roughpaper / quick tests
add_subdirectory(roughpaper EXCLUDE_FROM_ALL)
