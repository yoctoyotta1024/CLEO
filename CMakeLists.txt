cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) # if using Kokkos with NVC++ compiler 

project(project
  LANGUAGES CXX
  DESCRIPTION "CLEO by Clara Bayley and other developers"
  VERSION 4.0
)

add_subdirectory(libs)
add_subdirectory(src)