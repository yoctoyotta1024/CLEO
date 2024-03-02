# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("runCLEO")
message("${PROJECT_NAME} PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS gridboxes initialise observers runcleo superdrops zarr)

# create primary executable "runCLEO"
set(runCLEO "runCLEO")
add_executable(${runCLEO} "main.cpp")

# Add directories and libraries for target
target_link_libraries(${runCLEO} PRIVATE coupldyn_fromfile cartesiandomain ${SDMLIBS})
target_link_libraries(${runCLEO} PUBLIC Kokkos::kokkos)
target_include_directories(${runCLEO} PRIVATE "${CMAKE_SOURCE_DIR}/libs")
#target_include_directories(${runCLEO} PRIVATE ${PROJECT_SOURCE_DIR})

# set specific C++ compiler options for target (optional)
#target_compile_options(${runCLEO} PRIVATE)

# set C++ properties for target
set_target_properties(${runCLEO} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
