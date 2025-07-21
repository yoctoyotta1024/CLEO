# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("fromfile_irreg")
message(STATUS "CLEO including ${PROJECT_NAME} with PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# require MPI explicitly for this library
find_package(MPI REQUIRED COMPONENTS C)

# Set libraries from CLEO to link with executable
set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)

# create primary executable for CLEO in 2-D coupled to thermodynamics from file setup
add_executable(fromfile_irreg EXCLUDE_FROM_ALL "main_fromfile_irreg.cpp")

# Add directories and link libraries to target
target_include_directories(fromfile_irreg PRIVATE "${CLEO_SOURCE_DIR}/libs" ${MPI_INCLUDE_PATH})
target_link_libraries(fromfile_irreg PRIVATE coupldyn_fromfile cartesiandomain "${CLEOLIBS}")
target_link_libraries(fromfile_irreg PUBLIC Kokkos::kokkos MPI::MPI_C)

# set specific C++ compiler options for target (optional)
#target_compile_options(fromfile_irreg PRIVATE)

# set compiler properties for target(s)
set_target_properties(fromfile_irreg PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
