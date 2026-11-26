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


Running Examples on Different Machines:
----------------------------------------

Each example can be run by building Cleo, compiling the relevant executable, and then running the
example's Python script. There are bash scripts in ``scripts/`` to help you to do all this
relatively smoothly on DKRZ's Levante HPC, on JSC's JUPITER HPC, or on a generic/arbitrary,
so-called "vanilla", machine.


.. _how_bashscripts_work:

The Bash Scripts
#################

Each machine has its own bash scripts directory (``scripts/vanilla/``, ``scripts/levante/`` or
``scripts/jupiter/``). Scripts which are common to all the machines are kept in
``scripts/common/``.

The bash scripts for every machine are organised in the same way and your entry point is either
the ``cpu.sh`` or ``gpu.sh`` script inside ``scripts/[your_machine]/`` for your machine.
(Alternatively you can call ``build_compile_run_plot_cleo.sh`` directly,
see below: :ref:`build_compile_run_plot_cleo`.)

Inside ``scripts/[your_machine]/`` you will find:

* ``cpu.sh`` (and ``gpu.sh`` on an HPC): your entry point; the job script(s) you execute,
  or submit with ``sbatch``, to run one or more of the examples
* ``build_compile_run_plot_cleo.sh``: the script the job scripts (``cpu.sh`` and ``gpu.sh``)
  call for each example. It states which examples, build configurations and compilers your
  machine supports
* ``build_flags.sh`` and ``runtime_settings.sh``: the compiler flags, Kokkos flags and runtime
  environment for your machine
* ``helpers/``: the packages/modules used by your machine (and the helper script for you to
  :ref:`install YAC and YAXT<install_yac>`). It's highly unlikely that you need to familiarise
  yourself with scripts inside ``helpers/``

Additionally, every example is setup once, for all machines, in ``scripts/common/examples.sh``.
This is where you can find (or change) an example's default build directory, CMake flags,
executable(s), its Python script and the arguments given to that Python script.

Running An Example
####################

Running an example has five steps:

1) ``configure``: **configure Cleo's build** (with CMake in your chosen build directory)
2) ``compile``: **compile Cleo's executable(s)** (only the ones relevant to the chosen example)
3) ``inputs``: **generate the example's input files** (by calling the example's Python script)
4) ``run``: **run the executable(s)** (by calling the example's Python script)
5) ``plot``: **plot the results** (by calling the example's Python script)

You choose which steps a job script does with a comma-separated list, or ``all`` (the default)
for all five. The steps are always done in the order above. For example:

* ``all``: configure, compile, generate the input files, run and plot
* ``configure,compile``: e.g. on an HPC's login node before submitting a job to run the example
* ``compile,inputs,run,plot``: recompile, then generate the input files, run and plot, reusing an
  existing build
* ``compile,run``: recompile and run again without generating new input files
* ``plot``: only plot the results again

*Note*: the ``configure`` step needs internet access, because CMake downloads some of Cleo's
dependencies (e.g. Kokkos) when it configures Cleo. On machines whose compute nodes have no
internet access, such as JUPITER, first do ``configure,compile`` on a login node, and then submit a
job which does ``inputs,run,plot``.

To build an example from scratch, delete its build directory (e.g.
``rm -rf [CLEO_PATH2CLEO]/build_const2d``) and then do the ``configure`` step (e.g. with ``all``),
so that CMake configures Cleo again without reusing its old cache.


The Examples On Your Machine
###############################

Now you understand the bash scripts, you can try to run the examples which are
compatible with your machine:

.. toctree::
   :maxdepth: 1

   examples_vanilla
   examples_levante
   examples_jupiter


A little bit of extra information....
----------------------------------------

.. _build_compile_run_plot_cleo:

Using ``build_compile_run_plot_cleo.sh`` directly
##################################################

You can also skip the job scripts and call a machine's ``build_compile_run_plot_cleo.sh`` script
directly. It takes the same arguments on every machine:

.. code-block:: console

  $ scripts/[machine]/build_compile_run_plot_cleo.sh [example] [buildtype] [compilername] \
      [path2CLEO] [path2build] [build_flags] [yacyaxtroot] [enabledebug] \
      [stacksize_limit] [steps]

All the arguments are optional; use ``""`` to keep an argument's default. ``path2build`` is the
example's exact build directory (by default its directory inside your Cleo directory, e.g.
``build_colls0d/shima2009``) and
``steps`` is a comma-separated list of the steps above, or ``all``. E.g. to only configure and
compile the Shima et al. 2009 example with OpenMP and in debug mode on a vanilla machine, from
your Cleo directory:

.. code-block:: console

  $ scripts/vanilla/build_compile_run_plot_cleo.sh shima2009 openmp gcc $(pwd) /your/path/to/builds "" "" true "" configure,compile


The Bash Scripts in More Detail
#################################

You can learn more about the bash scripts :ref:`in Cleo's programming guide<bashscripts_detail>`.
