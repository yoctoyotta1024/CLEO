# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("bubble3d")
message(STATUS "CLEO including ${PROJECT_NAME} with PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)

# create primary executable for CLEO in 2-D coupled to thermodynamics from file setup
add_executable(bubble3d EXCLUDE_FROM_ALL "main_bubble3d.cpp")

# Add directories and link libraries to target
target_link_libraries(bubble3d PRIVATE coupldyn_yac cartesiandomain "${CLEOLIBS}")
target_link_libraries(bubble3d PUBLIC Kokkos::kokkos)
target_include_directories(bubble3d PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# set specific C++ compiler options for target (optional)
#target_compile_options(bubble3d PRIVATE)

# set compiler properties for target(s)
set_target_properties(bubble3d PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
