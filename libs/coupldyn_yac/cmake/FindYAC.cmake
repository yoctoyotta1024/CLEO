find_package(YAXT REQUIRED)
find_package(NetCDF REQUIRED)
find_package(LAPACK REQUIRED)
enable_language(C)
find_package(MPI REQUIRED COMPONENTS C)

if(YAXT_FOUND AND NetCDF_FOUND AND LAPACK_FOUND)
  # yac.h is YAC's supported C interface. Note yac_interface.h is the deprecated
  # interface and is only installed by YAC builds configured --enable-deprecated.
  find_path(YAC_C_INCLUDE_DIR
    NAMES yac.h
    DOC "YAC include dir")

  # YAC always installs its C libraries split into the message coupling interface
  # (libyac_mci.a) and the core (libyac_core.a). The combined libyac.a is only
  # produced by YAC builds configured --enable-deprecated. Link order matters:
  # libyac_mci.a depends on libyac_core.a.
  find_library(YAC_MCI_LIBRARY
    NAMES libyac_mci.a
    DOC "YAC C message coupling interface Library")

  find_library(YAC_CORE_LIBRARY
    NAMES libyac_core.a
    DOC "YAC C core Library")

  # bundled by YAC as libyac_mtime.a, or external (--with-external-mtime) as libmtime.a
  find_library(YAC_C_MTIME_LIBRARY
    NAMES libyac_mtime.a libmtime.a
    DOC "YAC C mtime Library")

  mark_as_advanced(YAC_C_INCLUDE_DIR
    YAC_MCI_LIBRARY
    YAC_CORE_LIBRARY
    YAC_C_MTIME_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(YAC
    REQUIRED_VARS YAC_MCI_LIBRARY YAC_CORE_LIBRARY YAC_C_INCLUDE_DIR YAC_C_MTIME_LIBRARY
  )

  if(YAC_FOUND)
    if(NOT TARGET YAC::YAC)
      add_library(YAC::YAC INTERFACE IMPORTED)
      target_include_directories(YAC::YAC INTERFACE "${YAC_C_INCLUDE_DIR}")
      target_link_libraries(YAC::YAC INTERFACE "${YAC_MCI_LIBRARY}" "${YAC_CORE_LIBRARY}" "${YAC_C_MTIME_LIBRARY}" YAXT::YAXT_C NetCDF::NetCDF_C MPI::MPI_C LAPACK::LAPACK fyaml m "-L${CLEO_FYAMLLIB}")
    endif()
  endif()
endif()
