.. _bashscripts_detail:

The Bash Scripts in Detail
==========================

This page is an optional read for anyone who wants to change the bash scripts in ``scripts/``,
e.g. to change the compiler flags, add a new example or add a new computer. To simply run the
examples, see :ref:`how the bash scripts work<bashscripts>` and the page for
:ref:`your computer<examples>` instead.

The scripts are designed so that everything which is the same on every computer is written once,
in ``scripts/common/``, and everything which is specific to a computer is in that computer's
directory. Every computer's directory has the same files, which provide the same functions, so
the common scripts can call them without needing to know which computer they are running on.


Layout
------

.. code-block:: text

  scripts/
  ├── common/                              the same on every computer
  │   ├── run_jobs.sh                      logic of the job scripts (modes, experiments list)
  │   ├── build_compile_run_plot_cleo.sh   the pipeline: checks, build, compile, run and plot
  │   ├── experiments.sh                   every example: build directory, CMake flags, executables,
  │   │                                    Python script and its arguments
  │   ├── build_cleo.sh                    configures Cleo with CMake
  │   ├── compile_cleo.sh                  compiles the executables with make
  │   ├── clean_build.sh                   deletes a build directory (make_clean)
  │   ├── build_openmp.sh                  Kokkos flags for OpenMP on the host
  │   ├── build_threads.sh                 Kokkos flags for C++ threads on the host
  │   ├── build_cuda.sh                    Kokkos flags for CUDA on the device
  │   ├── build_yac.sh                     CMake flags for YAC
  │   ├── mpi_compilers.sh                 chooses the mpic++ and mpicc compiler wrappers
  │   ├── check_inputs.sh                  input validation
  │   └── print_configuration.sh           prints the configuration of a run
  │
  └── [computer]/                          vanilla, levante or jupiter
      ├── cpu.sh (and gpu.sh)              the job script(s)
      ├── build_compile_run_plot_cleo.sh   what the computer supports, then calls the pipeline
      ├── build_flags.sh                   compiler flags and Kokkos flags
      ├── runtime_settings.sh              runtime environment
      └── helpers/
          ├── [computer]_packages.sh       modules/packages and their paths
          ├── build_yac.sh                 YAC flags (and compiler check) for this computer
          └── install_yac.sh               installs YAC and YAXT on this computer

Each script defines one or more functions which the other scripts call after sourcing it, e.g.
``scripts/common/experiments.sh`` defines ``load_experiment_config``. Many of them can also be
executed directly to call their function.


What Happens When You Run a Job Script
--------------------------------------

For example, ``sbatch scripts/levante/cpu.sh all constthermo2d`` does the following:

1) ``scripts/levante/cpu.sh`` sets the paths (``CLEO_PATH2CLEO``, ``CLEO_PYTHON``,
   ``CLEO_YACYAXTROOT`` and ``CLEO_PATH2BUILD``), the computer's environment (e.g. modules), its
   ``experiments`` list and ``CLEO_MAKE_CLEAN``, and then calls ``run_cleo_jobs`` from
   ``scripts/common/run_jobs.sh``.

2) ``run_cleo_jobs`` turns the mode into steps (``all``, ``build,compile`` or
   ``compile,run,plot``), checks the paths are not ``<...>`` placeholders, and then, for the given
   example or every entry of the ``experiments`` list, calls the computer's
   ``build_compile_run_plot_cleo.sh``.

3) ``scripts/levante/build_compile_run_plot_cleo.sh`` sets which build types, compilers and
   examples Levante supports and its defaults, and then calls ``build_compile_run_plot_cleo`` from
   ``scripts/common/build_compile_run_plot_cleo.sh``, which:

   a) reads the arguments, uses the computer's defaults for any which are empty, and checks them
      (and calls ``machine_check_inputs`` if the computer defines it),

   b) exports the ``CLEO_*`` environment variables and calls ``load_experiment_config`` from
      ``scripts/common/experiments.sh`` for the example's build directory, CMake flags,
      executables, Python script and arguments,

   c) prints the configuration,

   d) if ``make_clean`` is ``true``, deletes the example's build directory (``clean_cleo_build``),

   e) the ``build`` step: calls ``build_cleo`` from ``scripts/common/build_cleo.sh``, which calls
      ``configure_machine_build_flags`` (``levante/build_flags.sh``) and
      ``configure_machine_yac_flags`` (``levante/helpers/build_yac.sh``), and then configures Cleo
      with CMake,

   f) the ``compile`` step: calls ``compile_cleo`` from ``scripts/common/compile_cleo.sh``, which
      runs ``make`` for the example's executables. If there is no ``build`` step, it first calls
      ``configure_machine_build_flags`` to load the computer's compilers,

   g) the ``run`` and ``plot`` steps: call ``configure_machine_runtime_settings``
      (``levante/runtime_settings.sh``) and then the example's Python script, once with
      ``--do_inputfiles`` and once with ``--do_run_executable`` for the ``run`` step, and once
      with ``--do_plot_results`` for the ``plot`` step.

Inside ``configure_machine_build_flags``, the computer's modules are reset and loaded
(``[computer]_reset_modules`` and ``[computer]_load_build_stack`` from
``helpers/[computer]_packages.sh``), the compiler wrappers are chosen (``configure_mpi_compilers``),
the compiler flags and basic Kokkos flags are set, and the Kokkos flags for the build type are
added (``configure_openmp_build``, ``configure_threads_build`` or ``configure_cuda_build``).


.. _bashscripts_where:

Where Do I Change...?
---------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - To change...
     - Edit...
   * - the compiler flags (for debug and release builds)
     - ``configure_machine_build_flags`` in ``scripts/[computer]/build_flags.sh``, the
       ``CLEO_CXX_FLAGS`` for each compiler
   * - the Kokkos architecture and basic Kokkos flags
     - ``CLEO_KOKKOS_BASIC_FLAGS`` in ``scripts/[computer]/build_flags.sh``
   * - the Kokkos flags for a build type (e.g. CUDA)
     - ``scripts/common/build_openmp.sh``, ``build_threads.sh`` or ``build_cuda.sh`` for every
       computer, or the ``cuda)`` case in ``scripts/[computer]/build_flags.sh`` for one computer
   * - the modules or packages loaded, their versions and paths (e.g. to CUDA, NetCDF or fyaml)
     - ``scripts/[computer]/helpers/[computer]_packages.sh``
   * - the MPI compiler wrappers YAC must be built with
     - ``configure_machine_yac_flags`` in ``scripts/[computer]/helpers/build_yac.sh``
   * - the runtime environment (e.g. UCX, OpenMP thread placement, stack size)
     - ``configure_machine_runtime_settings`` in ``scripts/[computer]/runtime_settings.sh``
   * - which build types, compilers and examples a computer supports, its default build type,
       stack size and number of make jobs
     - the ``machine configuration`` section of
       ``scripts/[computer]/build_compile_run_plot_cleo.sh``
   * - an example's build directory, CMake flags, executables, Python script or the arguments to
       its Python script
     - the example's entry in ``load_experiment_config`` in ``scripts/common/experiments.sh``
   * - the CMake flags common to most examples
     - ``cleo_common_flags`` in ``scripts/common/experiments.sh``
   * - the Slurm options, which examples a job runs, the paths or ``CLEO_MAKE_CLEAN``
     - the job script, ``scripts/[computer]/cpu.sh`` or ``gpu.sh``
   * - the modes of the job scripts
     - ``run_cleo_jobs`` in ``scripts/common/run_jobs.sh``
   * - the steps, the order they run in or how the Python script is called
     - ``build_compile_run_plot_cleo`` in ``scripts/common/build_compile_run_plot_cleo.sh``

*Note*: if you change the compiler, build type or Kokkos flags of an existing build, build it from
scratch (``CLEO_MAKE_CLEAN=true``) so that CMake does not reuse its old cache.


Environment Variables
---------------------

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Variable
     - Set by
     - Meaning
   * - ``CLEO_MACHINE``
     - the job script and ``build_compile_run_plot_cleo.sh`` of each computer
     - the computer, i.e. its directory in ``scripts/``
   * - ``CLEO_PATH2CLEO``
     - the job script (the directory you run or submit it from), then the pipeline
     - path to your Cleo directory
   * - ``CLEO_PYTHON``
     - the job script (default ``[CLEO_PATH2CLEO]/.venv/bin/python3``)
     - Python interpreter which runs the examples' Python scripts
   * - ``CLEO_YACYAXTROOT``
     - the job script (a ``<...>`` placeholder you must replace)
     - directory containing your ``yac`` and ``yaxt`` installations
   * - ``CLEO_PATH2BUILD``
     - the job script (a ``<...>`` placeholder you must replace), then ``load_experiment_config``
     - in the job script, the directory in which examples are built; afterwards the example's own
       build directory, e.g. ``[CLEO_PATH2BUILD]/build_const2d/``
   * - ``CLEO_MAKE_CLEAN``
     - the job script (default ``false``)
     - ``true`` deletes each example's build directory before building it
   * - ``CLEO_MAKE_JOBS``
     - you (optional), else the computer's default
     - number of parallel ``make`` jobs
   * - ``CLEO_BUILDTYPE``, ``CLEO_COMPILERNAME``, ``CLEO_ENABLEDEBUG``
     - the pipeline, from its arguments
     - build type, compiler and whether to make a debug build
   * - ``CLEO_BUILD_FLAGS``
     - ``load_experiment_config``
     - the example's CMake flags
   * - ``CLEO_CXX_COMPILER``, ``CLEO_CC_COMPILER``
     - ``configure_mpi_compilers``
     - the MPI compiler wrappers
   * - ``CLEO_CXX_FLAGS``
     - ``configure_machine_build_flags``
     - the C++ compiler flags
   * - ``CLEO_KOKKOS_BASIC_FLAGS``, ``CLEO_KOKKOS_HOST_FLAGS``, ``CLEO_KOKKOS_DEVICE_FLAGS``
     - ``configure_machine_build_flags`` and the ``configure_[buildtype]_build`` functions
     - the Kokkos flags
   * - ``CLEO_CUDA_ROOT``
     - ``[computer]_load_build_stack`` for CUDA builds (or you)
     - path to the CUDA installation
   * - ``CLEO_YAC_FLAGS``
     - ``build_yac``
     - the CMake flags for YAC


Adding a New Example
--------------------

1. Add an entry for the example to ``load_experiment_config`` in ``scripts/common/experiments.sh``,
   setting its ``build_subdir``, ``build_flags``, ``executables``, ``pythonscript``,
   ``src_config_filename`` and ``script_args``. The Python script must accept the path to Cleo,
   the path to the build directory and ``script_args``, followed by one of ``--do_inputfiles``,
   ``--do_run_executable`` or ``--do_plot_results``, like the existing examples' Python scripts.

2. Add the example to ``machine_experiments`` in the ``build_compile_run_plot_cleo.sh`` of every
   computer which supports it.

3. Optionally, add it to the ``experiments`` list of a job script.


Adding a New Computer
---------------------

1. Copy the directory of the most similar computer, e.g. ``scripts/levante/`` for an HPC with Slurm
   and modules, or ``scripts/vanilla/`` for a computer without modules, and rename it and its
   ``helpers/[computer]_packages.sh`` file.

2. In ``helpers/[computer]_packages.sh``, set the computer's modules/packages and paths, and
   rename and adapt these functions (on a computer without modules they can do nothing, as for the
   vanilla computer):

   * ``[computer]_reset_modules``: unloads all modules,
   * ``[computer]_load_build_stack``: loads the compilers, MPI and CMake for a compiler and build
     type (and sets ``CLEO_CUDA_ROOT`` for CUDA),
   * ``[computer]_load_runtime_stack``: loads what the executables need to run,
   * ``[computer]_load_yac_dependencies``: loads what YAC needs,
   * ``[computer]_fyamllib_for_compiler``: prints the path to the fyaml library for a compiler.

3. In ``build_flags.sh``, ``runtime_settings.sh`` and ``helpers/build_yac.sh``, replace the
   old computer's name in the calls to these functions, and set the computer's compiler flags,
   Kokkos flags, runtime environment and YAC compiler check. Each of these files must still
   define ``configure_machine_build_flags``, ``configure_machine_runtime_settings`` and
   ``configure_machine_yac_flags`` respectively, because the common scripts call these functions.

4. In ``build_compile_run_plot_cleo.sh``, set ``CLEO_MACHINE`` to the new directory's name and set
   the ``machine configuration``: ``machine_default_buildtype``, ``machine_buildtypes``,
   ``machine_compilers``, ``machine_experiments``, ``machine_default_stacksize`` and
   ``machine_default_make_jobs``. Optionally define ``machine_check_inputs`` for any extra checks.

5. In the job script(s), set ``CLEO_MACHINE``, the Slurm options, the environment and the
   ``experiments`` list.

6. Update ``helpers/install_yac.sh`` for the computer, and add a page for it to these docs.


Safety Checks and Things to Know
--------------------------------

* The job scripts stop with an error if any of the paths are still ``<...>`` placeholders.

* ``CLEO_PATH2CLEO`` is the directory you run or submit a job script from, so always run or submit
  the job scripts from your Cleo directory.

* When building from scratch, only build directories which contain a ``CMakeCache.txt`` file are
  deleted, and never ``/``, your home directory, your Cleo directory or a directory containing it.

* ``make_clean=true`` needs the ``build`` step, and the ``run`` mode reuses an existing build, so
  building from scratch is not possible in the ``run`` mode.

* If you choose an example when running a job script but not its build type, the computer's
  default build type is used, also for ``gpu.sh``. So for GPUs also give the ``cuda`` build type.

* ``CLEO_YACYAXTROOT`` must be set even if no example needs YAC.
