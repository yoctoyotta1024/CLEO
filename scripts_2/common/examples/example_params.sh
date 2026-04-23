#!/bin/bash

### Central experiment parameter lookup.
### Set 'experiment' before sourcing this file, or pass the experiment name as $3.
### Positional args allow CLI overrides:
###   $1 = path2build  override (optional, leave empty to use default)
###   $2 = build_flags override (optional, leave empty to use default)
###   $3 = experiment  name     (optional, overrides the 'experiment' variable)
###

path2build_override=${1:-""}
build_flags_override=${2:-""}
experiment=${3:-${experiment}}

if [ "${path2CLEO}" == "" ]; then
  echo "Please provide path to CLEO source directory"
  exit 1
fi

if [[ -z "${experiment}" ]]; then
  echo "Error: 'experiment' must be set before sourcing experiment_params.sh (or passed as \$4)"
  exit 1
fi

case "${experiment}" in

  as2017)
    path2build=${path2CLEO}/build_adia0d/as2017
    build_flags="-DCLEO_COUPLED_DYNAMICS=cvode -DCLEO_DOMAIN=cartesian \
      -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
    executables="adia0d"

    pythonscript=${path2CLEO}/examples/adiabaticparcel/as2017.py
    src_config_filename=${path2CLEO}/examples/adiabaticparcel/src/config/as2017_config.yaml
    script_args="${src_config_filename} \
      --do_inputfiles --do_run_executable --do_plot_results"

    ;;

  breakup)
      path2build=${path2CLEO}/build_colls0d/breakup/
      build_flags="-DCLEO_COUPLED_DYNAMICS=null -DCLEO_DOMAIN=cartesian \
        -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
      executables="longcolls lowlistcolls szakallurbichcolls testikstraubcolls"

      pythonscript=${path2CLEO}/examples/boxmodelcollisions/breakup.py
      src_config_filename=${path2CLEO}/examples/boxmodelcollisions/src/config/breakup_config.yaml
      script_args="${src_config_filename} --kernels long lowlist szakallurbich testikstraub \
        --do_inputfiles --do_run_executable --do_plot_results"

      ;;

  constthermo2d)
      path2build=${path2CLEO}/build_const2d/
      build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
        -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
      executables="const2d"

      pythonscript=${path2CLEO}/examples/constthermo2d/constthermo2d.py
      src_config_filename=${path2CLEO}/examples/constthermo2d/src/config/const2d_config.yaml
      script_args="${src_config_filename} \
        --do_inputfiles --do_run_executable --do_plot_results"
      ;;

  cuspbifurc)
    path2build=${path2CLEO}/build_adia0d/cuspbifurc/
    build_flags="-DCLEO_COUPLED_DYNAMICS=cvode -DCLEO_DOMAIN=cartesian \
      -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
    executables="adia0d"

    pythonscript=${path2CLEO}/examples/adiabaticparcel/cuspbifurc.py
    src_config_filename=${path2CLEO}/examples/adiabaticparcel/src/config/cuspbifurc_config.yaml
    script_args="${src_config_filename} \
      --do_inputfiles --do_run_executable --do_plot_results"
    ;;

  divfree2d)
    path2build=${path2CLEO}/build_divfree2d/
    build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
      -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
    executables="divfree2d"

    pythonscript=${path2CLEO}/examples/divfreemotion/divfree2d.py
    src_config_filename=${path2CLEO}/examples/divfreemotion/src/config/divfree2d_config.yaml
    script_args="${src_config_filename} \
      --do_inputfiles --do_run_executable --do_plot_results"
    ;;

  eurec4a1d)
    path2build=${path2CLEO}/build_eurec4a1d/
    build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
      -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
    executables="eurec4a1d"

    pythonscript=${path2CLEO}/examples/eurec4a1d/eurec4a1d.py
    src_config_filename=${path2CLEO}/examples/eurec4a1d/src/config/eurec4a1d_config.yaml
    script_args="${src_config_filename} \
      --do_inputfiles --do_run_executable --do_plot_results"
    ;;

  python_bindings)
  path2build=${path2CLEO}/build_pybind/
  build_flags="-DCLEO_COUPLED_DYNAMICS=numpy -DCLEO_DOMAIN=cartesian \
    -DCLEO_NO_ROUGHPAPER=true -DCLEO_PYTHON=${CLEO_PYTHON}"
  executables="cleo_python_bindings"

  pythonscript=${path2CLEO}/examples/python_bindings/python_bindings.py
  src_config_filename=${path2CLEO}/examples/python_bindings/src/config/pybind_config.yaml
  script_args="${src_config_filename} \
    --do_inputfiles --do_run_executable --do_plot_results"
  ;;

  rainshaft1d)
    path2build=${path2CLEO}/build_rshaft1d/
    build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
      -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
    executables="rshaft1d"

    pythonscript=${path2CLEO}/examples/rainshaft1d/rainshaft1d.py
    src_config_filename=${path2CLEO}/examples/rainshaft1d/src/config/rshaft1d_config.yaml
    script_args="${src_config_filename} \
      --do_inputfiles --do_run_executable --do_plot_results"
    ;;

  shima2009)
    path2build=${path2CLEO}/build_colls0d/shima2009/
    build_flags="-DCLEO_COUPLED_DYNAMICS=null -DCLEO_DOMAIN=cartesian \
      -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
    executables="golcolls longcolls"

    pythonscript=${path2CLEO}/examples/boxmodelcollisions/shima2009.py
    src_config_filename=${path2CLEO}/examples/boxmodelcollisions/src/config/shima2009_config.yaml
    script_args="${src_config_filename} --kernels golovin long1 long2 \
      --do_inputfiles --do_run_executable --do_plot_results"
    ;;
  *)
    echo "Error: unknown experiment '${experiment}'."
    echo "Available: as2017 cuspbifurc breakup shima2009 constthermo2d divfree2d"
    echo "           eurec4a1d rainshaft1d python_bindings kokkostools"
    echo "           fromfile fromfile_irreg bubble3d"
    exit 1
    ;;
esac

### Apply CLI overrides if provided (non-empty args take precedence over defaults)
[[ -n "${path2build_override}" ]]  && path2build="${path2build_override}"
[[ -n "${build_flags_override}" ]] && build_flags="${build_flags_override}"

export CLEO_PATH2BUILD=${path2build}
export CLEO_BUILD_FLAGS="${build_flags}"
