# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("roughpaper_src_cleocoupledsdm")
message(STATUS "CLEO including ${PROJECT_NAME} with PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS configuration gridboxes initialise observers runcleo superdrops zarr)

# create primary executable "cleocoupledsdm" for cleo sdm coupled to dynamics solver
add_executable(cleocoupledsdm main.cpp)

# Add directories and link libraries to target
target_link_libraries(cleocoupledsdm PRIVATE coupldyn_fromfile cartesiandomain ${CLEOLIBS})
target_link_libraries(cleocoupledsdm PUBLIC Kokkos::kokkos)
target_include_directories(cleocoupledsdm PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory
target_include_directories(cleocoupledsdm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# set specific C++ compiler options for target (optional)
#target_compile_options(cleocoupledsdm PRIVATE)

# set C++ properties for target
set_target_properties(cleocoupledsdm PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
