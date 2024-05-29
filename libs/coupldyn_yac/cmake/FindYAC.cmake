find_package(YAXT)
find_package(NetCDF)
find_package(LAPACK)
enable_language(C)
find_package(MPI REQUIRED COMPONENTS C)

if(YAXT_FOUND AND NetCDF_FOUND AND LAPACK_FOUND)
  find_path(YAC_C_INCLUDE_DIR
    NAMES yac_interface.h
    DOC "YAC include dir")

  find_path(YAC_FORTRAN_INCLUDE_DIR
    NAMES mo_yac_finterface.mod
    HINTS ${YAC_C_INCLUDE_DIR}
    DOC "YAC fortran include dir")

  find_library(YAC_C_LIBRARY
    NAMES libyac.a
    DOC "YAC C Library")

  find_library(YAC_C_MTIME_LIBRARY
    NAMES libyac_mtime.a libmtime.a
    DOC "YAC C mtime Library")

  mark_as_advanced(YAC_C_INCLUDE_DIR
    YAC_C_LIBRARY
    YAC_C_MTIME_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(YAC
    REQUIRED_VARS YAC_C_LIBRARY YAC_C_INCLUDE_DIR YAC_C_MTIME_LIBRARY
  )

  if(YAC_FOUND)
    if(NOT TARGET YAC::YAC)
      add_library(YAC::YAC INTERFACE IMPORTED)
      target_include_directories(YAC::YAC INTERFACE "${YAC_C_INCLUDE_DIR}" "${YAC_FORTRAN_INCLUDE_DIR}")
      target_link_libraries(YAC::YAC INTERFACE "${YAC_C_LIBRARY}" "${YAC_C_MTIME_LIBRARY}" YAXT::YAXT_C NetCDF::NetCDF_C MPI::MPI_C LAPACK::LAPACK fyaml m "-L/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib")
    endif()
  endif()
endif()
