# Finds YAC through its pkg-config file yac-mci.pc, which requires yac-core.pc.
# Together they list YAC's libraries in link order, including the transitive
# dependencies (NetCDF, YAXT, LAPACK, mtime, fyaml). YAC_ROOT is searched first.
#
# Defines the imported target YAC::YAC and the variable YAC_C_INCLUDE_DIR.
#
# MPI is linked separately because the .pc files omit it when YAC was built
# with CC=mpicc. A .pc file whose prefix no longer matches the install location
# (e.g. a relocated install) can be used by passing
# -DPKG_CONFIG_ARGN=--define-prefix (CMake >= 3.22).

enable_language(C)
find_package(MPI REQUIRED COMPONENTS C)
find_package(PkgConfig REQUIRED)

set(_yac_saved_prefix_path "${CMAKE_PREFIX_PATH}")
if(YAC_ROOT)
  list(PREPEND CMAKE_PREFIX_PATH "${YAC_ROOT}")
endif()
pkg_check_modules(YAC_PC QUIET IMPORTED_TARGET yac-mci)
set(CMAKE_PREFIX_PATH "${_yac_saved_prefix_path}")
unset(_yac_saved_prefix_path)

if(YAC_PC_FOUND)
  # yac.h is YAC's supported C interface (yac_interface.h is deprecated)
  find_path(YAC_C_INCLUDE_DIR
    NAMES yac.h
    PATHS "${YAC_PC_INCLUDEDIR}"
    NO_DEFAULT_PATH
    DOC "YAC include dir")
  mark_as_advanced(YAC_C_INCLUDE_DIR)

  # YAC's own CMake build writes incomplete .pc files in these versions
  if(NOT YAC_PC_LIBRARIES MATCHES "mtime")
    message(FATAL_ERROR "yac-mci.pc of YAC ${YAC_PC_VERSION} does not link mtime. "
      "This is a bug of YAC 3.20.x built with CMake; use YAC >= 3.21 or build YAC with autotools.")
  endif()
  if(YAC_PC_VERSION VERSION_LESS 3.19 AND EXISTS "${YAC_PC_LIBDIR}/cmake/yac/yac-config.cmake")
    message(FATAL_ERROR "yac-core.pc of YAC ${YAC_PC_VERSION} built with CMake does not link LAPACK. "
      "This is a bug of YAC 3.15 to 3.18 built with CMake; use YAC >= 3.19 or build YAC with autotools.")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAC
  REQUIRED_VARS YAC_C_INCLUDE_DIR YAC_PC_FOUND
  VERSION_VAR YAC_PC_VERSION
)

if(YAC_FOUND AND NOT TARGET YAC::YAC)
  add_library(YAC::YAC INTERFACE IMPORTED)
  target_link_libraries(YAC::YAC INTERFACE PkgConfig::YAC_PC MPI::MPI_C)
endif()
