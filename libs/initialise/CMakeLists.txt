# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library)
set(LIBNAME "initialise")
message(STATUS "CLEO using ${LIBNAME} with LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIB_BINARY_DIR ${CMAKE_BINARY_DIR}/lib)

# require MPI explicitly for this library
find_package(MPI REQUIRED COMPONENTS C)

# Add executables and create library target
set(SOURCES
"gbx_bounds_from_binary.cpp"
"init_all_supers_from_binary.cpp"
"init_supers_from_binary.cpp"
"readbinary.cpp"
"timesteps.cpp"
)
# must use STATIC (not(!) SHARED) lib for linking to executable if build is CUDA enabled with Kokkos
add_library("${LIBNAME}" STATIC ${SOURCES})

# Add directories for target library
#target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PRIVATE "${CLEO_SOURCE_DIR}/libs" ${MPI_INCLUDE_PATH})

# Link libraries to target library
set(LINKLIBS configuration)
target_link_libraries(${LIBNAME} PUBLIC ${LINKLIBS})
target_link_libraries(${LIBNAME} PUBLIC Kokkos::kokkos MPI::MPI_C)

# set specific C++ compiler options for target (optional)
#target_compile_options(${LIBNAME} PRIVATE)

# set C++ properties for target
set_target_properties(${LIBNAME} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  ARCHIVE_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
