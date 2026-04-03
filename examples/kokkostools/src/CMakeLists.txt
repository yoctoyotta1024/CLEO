# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("kokkostools")
message(STATUS "CLEO including ${PROJECT_NAME} with PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)

# create primary executable for CLEO in serial, gpu and cpu kokkostools example test executable
add_executable(kokkostools EXCLUDE_FROM_ALL "main_kokkostools.cpp")

# Add directories and link libraries to target
target_link_libraries(kokkostools PRIVATE coupldyn_fromfile cartesiandomain "${CLEOLIBS}")
target_link_libraries(kokkostools PUBLIC Kokkos::kokkos)
target_include_directories(kokkostools PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# set specific C++ compiler options for target (optional)
#target_compile_options(kokkostools PRIVATE)

# set compiler properties for target(s)
set_target_properties(kokkostools PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
