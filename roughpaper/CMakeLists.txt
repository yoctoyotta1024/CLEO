# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("test")
message("${PROJECT_NAME} PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Set libraries from CLEO to link with executable
set(SDMLIBS superdrops)

# create primary executable "test"
set(test "test")
add_executable(${test} "main.cpp")

# Add directories and link libraries to target
target_link_libraries(${test} PRIVATE "${SDMLIBS}")
target_link_libraries(${test} PUBLIC Kokkos::kokkos)
target_compile_options(${test} PRIVATE ${compileops})
target_include_directories(${test} PRIVATE "${CMAKE_SOURCE_DIR}/libs")

# set specific C++ compiler options for target (optional)
target_compile_options(${test} PRIVATE -Werror -Wall -pedantic -g -gdwarf-4 -O0 -mpc64)

# set C++ properties for target
set_target_properties(${runCLEO} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
