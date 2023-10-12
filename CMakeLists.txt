# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# e.g. could run cmake with these compiler options:
# CXX=g++ CC=gcc cmake -S ./ -B ./build -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON

# set cmake version
cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) ### if using Kokkos with NVC++ compiler on Levante

project(project
  LANGUAGES CXX
  DESCRIPTION "CLEO by Clara Bayley and other developers"
  VERSION 1.0
)

add_subdirectory(libs)
add_subdirectory(src)