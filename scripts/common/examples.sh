#!/bin/bash

### Central example parameter lookup.
### Usage: load_example_config <example> [path2build] [build_flags]
###   example      name of the example
###   path2build   exact build directory (optional, empty: <CLEO_PATH2CLEO>/build_xxx)
###   build_flags  CMake flags replacing the example's own (optional)
### Requires CLEO_PATH2CLEO to be exported. Exports CLEO_PATH2BUILD and
### CLEO_BUILD_FLAGS, and sets executables, pythonscript and script_args.

load_example_config() {
  example=$1
  local path2build_override=${2:-}
  local build_flags_override=${3:-}

  ### Common CMake flags shared by most examples
  cleo_common_flags="-DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"

  ### build_subdir below is the example's default build directory inside CLEO_PATH2CLEO,
  ### used only when no exact build directory (path2build_override) is given.

  case "${example}" in

    as2017)
      build_subdir=build_adia0d/as2017
      build_flags="-DCLEO_COUPLED_DYNAMICS=cvode -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="adia0d"

      pythonscript=${CLEO_PATH2CLEO}/examples/adiabaticparcel/as2017.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/adiabaticparcel/src/config/as2017_config.yaml
      script_args="${src_config_filename}"
      ;;

    breakup)
      build_subdir=build_colls0d/breakup/
      build_flags="-DCLEO_COUPLED_DYNAMICS=null -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="longcolls lowlistcolls szakallurbichcolls testikstraubcolls"

      pythonscript=${CLEO_PATH2CLEO}/examples/boxmodelcollisions/breakup.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/boxmodelcollisions/src/config/breakup_config.yaml
      script_args="${src_config_filename} --kernels long lowlist szakallurbich testikstraub"
      ;;

    bubble3d)
      build_subdir=build_bubble3d/
      build_flags="-DCLEO_COUPLED_DYNAMICS=yac -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="bubble3d"
      path2experiment="/work/mh0731/m300950/icon/icon/experiments/bubble_cleo"
      path2iconfiles="/work/mh0731/m300950/icon/icon/experiments/bubble_1mom/outdata"

      pythonscript=${CLEO_PATH2CLEO}/examples/bubble3d/bubble3d.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml
      src_yac_config_filename=${CLEO_PATH2CLEO}/examples/bubble3d/src/config/yac_icon_cleo_coupling_config.yaml
      script_args="${src_config_filename} ${src_yac_config_filename} ${path2experiment} ${path2iconfiles}"
      ;;

    constthermo2d)
      build_subdir=build_const2d/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="const2d"

      pythonscript=${CLEO_PATH2CLEO}/examples/constthermo2d/constthermo2d.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/constthermo2d/src/config/const2d_config.yaml
      script_args="${src_config_filename}"
      ;;

    cuspbifurc)
      build_subdir=build_adia0d/cuspbifurc/
      build_flags="-DCLEO_COUPLED_DYNAMICS=cvode -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="adia0d"

      pythonscript=${CLEO_PATH2CLEO}/examples/adiabaticparcel/cuspbifurc.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/adiabaticparcel/src/config/cuspbifurc_config.yaml
      script_args="${src_config_filename}"
      ;;

    divfree2d)
      build_subdir=build_divfree2d/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="divfree2d"

      pythonscript=${CLEO_PATH2CLEO}/examples/divfreemotion/divfree2d.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/divfreemotion/src/config/divfree2d_config.yaml
      script_args="${src_config_filename}"
      ;;

    eurec4a1d)
      build_subdir=build_eurec4a1d/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="eurec4a1d"

      pythonscript=${CLEO_PATH2CLEO}/examples/eurec4a1d/eurec4a1d.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/eurec4a1d/src/config/eurec4a1d_config.yaml
      script_args="${src_config_filename}"
      ;;

    fromfile)
      build_subdir=build_fromfile/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="fromfile"

      pythonscript=${CLEO_PATH2CLEO}/examples/fromfile/fromfile.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/fromfile/src/config/fromfile_config.yaml
      script_args="${src_config_filename} --ntasks=4"
      ;;

    fromfile_irreg)
      build_subdir=build_fromfile_irreg/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="fromfile_irreg"

      pythonscript=${CLEO_PATH2CLEO}/examples/fromfile_irreg/fromfile_irreg.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/fromfile_irreg/src/config/fromfile_irreg_config.yaml
      script_args="${src_config_filename} --ntasks=4"
      ;;

    python_bindings)
      build_subdir=build_pybind/
      # pybind11 turns on link-time optimisation (LTO) for Release builds unless
      # CMAKE_INTERPROCEDURAL_OPTIMIZATION is set; LTO fails to link CUDA (nvcc) objects
      build_flags="-DCLEO_COUPLED_DYNAMICS=numpy -DCLEO_DOMAIN=cartesian \
        -DCLEO_NO_ROUGHPAPER=true -DCLEO_PYTHON=${CLEO_PYTHON} \
        -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF"
      executables="cleo_python_bindings"

      pythonscript=${CLEO_PATH2CLEO}/examples/python_bindings/python_bindings.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/python_bindings/src/config/pybind_config.yaml
      script_args="${src_config_filename}"
      ;;

    rainshaft1d)
      build_subdir=build_rshaft1d/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="rshaft1d"

      pythonscript=${CLEO_PATH2CLEO}/examples/rainshaft1d/rainshaft1d.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/rainshaft1d/src/config/rshaft1d_config.yaml
      script_args="${src_config_filename}"
      ;;

    roughpaper)
      # your own executable: roughpaper/src/main.cpp (note: roughpaper must not be
      # excluded, so cleo_common_flags with -DCLEO_NO_ROUGHPAPER=true is not used)
      build_subdir=build_roughpaper/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian -DCLEO_NO_PYBINDINGS=true"
      executables="cleocoupledsdm"

      pythonscript=${CLEO_PATH2CLEO}/roughpaper/roughpaper.py
      src_config_filename=${CLEO_PATH2CLEO}/roughpaper/src/config/config.yaml
      script_args="${src_config_filename}"
      ;;

    shima2009)
      build_subdir=build_colls0d/shima2009/
      build_flags="-DCLEO_COUPLED_DYNAMICS=null -DCLEO_DOMAIN=cartesian ${cleo_common_flags}"
      executables="golcolls longcolls"

      pythonscript=${CLEO_PATH2CLEO}/examples/boxmodelcollisions/shima2009.py
      src_config_filename=${CLEO_PATH2CLEO}/examples/boxmodelcollisions/src/config/shima2009_config.yaml
      script_args="${src_config_filename} --kernels golovin long1 long2"
      ;;

    *)
      echo "Error: unknown example '${example}'."
      exit 1
      ;;
  esac

  path2build=${path2build_override:-${CLEO_PATH2CLEO%/}/${build_subdir}}

  ### Apply CLI overrides if provided (non-empty args take precedence over defaults)
  [[ -n "${build_flags_override}" ]] && build_flags="${build_flags_override}"

  export CLEO_PATH2BUILD="${path2build}"
  export CLEO_BUILD_FLAGS="${build_flags}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  load_example_config "$@"
fi
