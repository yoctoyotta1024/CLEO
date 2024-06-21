# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("adiabaticparcel")
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS gridboxes initialise observers runcleo superdrops zarr)

# create primary executable for CVODE-CLEO coupled model adiabatic parcel setup
add_executable(adia0D EXCLUDE_FROM_ALL "main_adia0D.cpp")

# Add directories and link libraries to target
target_link_libraries(adia0D PRIVATE coupldyn_cvode cartesiandomain "${CLEOLIBS}")
target_link_libraries(adia0D PUBLIC Kokkos::kokkos)
target_include_directories(adia0D PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# set specific C++ compiler options for target (optional)
#target_compile_options(adia0D PRIVATE)

# set compiler properties for target(s)
set_target_properties(adia0D PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
