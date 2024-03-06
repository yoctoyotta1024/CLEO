# module load gcc && CXX=g++ CC=gcc cmake -S ./ -B ./build -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON

# set cmake version
cmake_minimum_required(VERSION 3.21.1)

project(project
  LANGUAGES CXX
  DESCRIPTION "roughpaper test space"
  VERSION 1.0
)

message("${PROJECT_NAME} PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

#set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic")

# add kokkos
add_subdirectory(extern/kokkos)

# add local sources
add_subdirectory(src)

# create primary executable "test"
set(test "test")
add_executable(${test} "main.cpp")

# Add directories and link libraries to target
target_link_libraries(${test} PUBLIC "${SDMLIBS}")
target_link_libraries(${test} PUBLIC Kokkos::kokkos)
target_include_directories(${test} PUBLIC "${CMAKE_SOURCE_DIR}/src")

# set C++ properties for target
set_target_properties(${test} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
