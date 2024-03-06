# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.21.1)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library)
message("Library 'src' LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# Add executables and create library target
add_library(src collisionkinetics.cpp)

# Add directories and link libraries for target library
target_include_directories(src PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(src PUBLIC Kokkos::kokkos)

# set C++ properties for target
set_target_properties(src PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
