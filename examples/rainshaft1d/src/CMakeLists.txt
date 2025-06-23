# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("rainshaft1d")
message(STATUS "CLEO including ${PROJECT_NAME} with PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)

# create primary executable for CLEO in 1-D coupled to thermodynamics from file setup
add_executable(rshaft1d EXCLUDE_FROM_ALL "main_rshaft1d.cpp")

# Add directories and link libraries to target
target_link_libraries(rshaft1d PRIVATE coupldyn_fromfile cartesiandomain "${CLEOLIBS}")
target_link_libraries(rshaft1d PUBLIC Kokkos::kokkos)
target_include_directories(rshaft1d PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# set specific C++ compiler options for target (optional)
#target_compile_options(rshaft1d PRIVATE)

# set compiler properties for target(s)
set_target_properties(rshaft1d PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
