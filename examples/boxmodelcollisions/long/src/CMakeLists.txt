# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  # cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("boxmodelcollisions_long")
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS gridboxes initialise observers runcleo superdrops zarr)

# create primary executable for CLEO uncoupled model long kernel box model setup
add_executable(longcolls EXCLUDE_FROM_ALL "main_longcolls.cpp")

# Add directories and link libraries to target
target_link_libraries(longcolls PRIVATE cartesiandomain "${CLEOLIBS}")
target_link_libraries(longcolls PUBLIC Kokkos::kokkos)
target_include_directories(longcolls PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory

# set specific C++ compiler options for target (optional)
#target_compile_options(longcolls PRIVATE)

# set compiler properties for target(s)
set_target_properties(longcolls PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
