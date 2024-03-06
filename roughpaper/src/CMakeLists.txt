# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.21.1)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library)
set(LIBNAME "src")
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR})

# Add executables and create library target
set(SOURCES
  "collisionkinetics.cpp"
)
add_library("${LIBNAME}" SHARED ${SOURCES})

# Add directories and link libraries for target library
target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${LIBNAME} PUBLIC Kokkos::kokkos)

# set C++ properties for target
set_target_properties(${LIBNAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
