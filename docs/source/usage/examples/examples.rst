.. _examples:

Examples
========

There are various examples of Cleo, with different build configurations, domains, microphysics,
coupling, and super-droplet motion etc. They can be found in the ``CLEO/examples`` directory. If you
would like to a copy of the reference solutions please :ref:`contact us <contact>`.


.. _install_plotcleo:

Installing ``plotcleo``:
------------------------

Before being able to run the examples you will need to locally install the ``plotcleo`` python
package from the ``examples/exampleplotting/`` directory. E.g.

.. code-block:: console

  $ uv build examples/exampleplotting/plotcleo
  $ uv pip install examples/exampleplotting/plotcleo/dist/plotcleo-[version].tar.gz
  $ uv run python -c "import plotcleo"


Running Examples on Different Computers:
----------------------------------------

Each example can be run by building Cleo, compiling the relevant executable, and then running the
example's Python script. There are bash scripts in ``scripts_2/`` to help you to do all this
relatively smoothly on DKRZ's Levante HPC, on JSC's JUPITER HPC, or on a generic/arbitrary,
so-called "vanilla", computer:

.. toctree::
   :maxdepth: 1

   examples_vanilla
   examples_levante
   examples_jupiter


.. _bashscripts:

How the Bash Scripts Work:
--------------------------

The bash scripts are organised in the same way for every computer. The steps which are the same
on all computers are in ``scripts_2/common/``, and each computer has its own directory
(``scripts_2/vanilla/``, ``scripts_2/levante/`` or ``scripts_2/jupiter/``) containing:

* ``cpu.sh`` (and ``gpu.sh`` on an HPC): the job script(s) you execute, or submit with ``sbatch``,
  to run one or more of the examples,
* ``build_compile_run_plot_cleo.sh``: the script the job scripts call for each example. It states
  which examples, build configurations and compilers that computer supports,
* ``build_flags.sh`` and ``runtime_settings.sh``: the compiler flags, Kokkos flags and runtime
  environment for that computer,
* ``helpers/``: the packages (modules) used on that computer, and a script to
  :ref:`install YAC and YAXT<install_yac>`.

Every example is described once, for all computers, in ``scripts_2/common/experiments.sh``. This is
where you can find (or change) an example's build directory, CMake flags, executable(s),
Python script and the arguments given to its Python script.

Running an example has four steps:

1) ``build``: configure Cleo with CMake in the example's build directory,
2) ``compile``: compile the example's executable(s),
3) ``run``: generate the input files and run the executable(s) by calling the example's Python script,
4) ``plot``: plot the results by calling the example's Python script again.

The job scripts combine these steps into three modes:

* ``all`` (the default): all four steps,
* ``build``: only the ``build`` and ``compile`` steps, e.g. on an HPC's login node before submitting
  a job to run the example,
* ``run``: recompile, then run and plot, reusing an existing build.

To build an example from scratch, i.e. to delete its build directory before configuring Cleo
with CMake again, set ``CLEO_MAKE_CLEAN=true`` when using the ``all`` or ``build`` modes. For
safety, only build directories which contain a ``CMakeCache.txt`` file are deleted.

.. dropdown:: Using ``build_compile_run_plot_cleo.sh`` directly
  :animate: fade-in

  You can also skip the job scripts and call a computer's ``build_compile_run_plot_cleo.sh`` script
  directly. It takes the same arguments on every computer:

  .. code-block:: console

    $ scripts_2/[computer]/build_compile_run_plot_cleo.sh [experiment] [buildtype] [compilername] \
        [path2CLEO] [path2build] [build_flags] [yacyaxtroot] [enabledebug] [make_clean] \
        [stacksize_limit] [steps]

  All the arguments are optional; use ``""`` to keep an argument's default. ``path2build`` is the
  directory in which each example's build directory is made (by default your Cleo directory) and
  ``steps`` is a comma-separated list of the steps above, or ``all``. E.g. to only build and
  compile the Shima et al. 2009 example with OpenMP and in debug mode on a vanilla computer, from
  your Cleo directory:

  .. code-block:: console

    $ scripts_2/vanilla/build_compile_run_plot_cleo.sh shima2009 openmp gcc $(pwd) /your/path/to/builds "" "" true "" "" build,compile
