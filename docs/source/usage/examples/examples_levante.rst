.. _examples_levante:

Examples on Levante
===================

Having :ref:`installed plotcleo<install_plotcleo>`, the following instructions are intended to guide you
through running each example using the bash scripts in ``scripts_2/levante/``. See
:ref:`how the bash scripts work<bashscripts>` for an overview of these scripts.

*Note*: the GPU job script, ``scripts_2/levante/gpu.sh``, chooses a build configuration which uses
GPUs. It must therefore run on a node in the GPU partition of Levante
(`see here <https://docs.dkrz.de/doc/levante/running-jobs/partitions-and-limits.html>`_
for documentation on Levante's partitions).

.. _configurebash_levante:

Configure the Bash Scripts
--------------------------

Every example is run with the same job script, ``scripts_2/levante/cpu.sh``
(or ``scripts_2/levante/gpu.sh`` to use GPUs).
Before using it for the first time, you will need to set the paths in the section at the top of
the script marked ``paths (EDIT THESE FOR YOUR SITE)``:

.. code-block:: console

  export CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"
  export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"
  export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-<PATH/TO/YACYAXT/INSTALL>}"
  export CLEO_PATH2BUILD="${CLEO_PATH2BUILD:-<PATH/TO/BUILD/ROOT>}"

You will need to configure ``cpu.sh`` and ``gpu.sh`` in the following ways:

* Set the path to your YAC and YAXT installations:

  replace ``<PATH/TO/YACYAXT/INSTALL>`` with the path to the directory containing your yac and yaxt
  directories (see :ref:`how to install YAC and YAXT<install_yac>`). If you do not intend to run an
  example that requires YAC, the path is not used, but it must still be set.

* Choose your build directory:

  replace ``<PATH/TO/BUILD/ROOT>`` with the directory in which you want Cleo to be built. Each
  example is built in its own directory inside this one, e.g. the Arabas and Shima 2017 example is
  built in ``<PATH/TO/BUILD/ROOT>/build_adia0d/as2017/``. (*hint*: to build in your Cleo directory
  use ``${CLEO_PATH2CLEO}``.)

* Use your Python version:

  by default ``CLEO_PYTHON`` is the Python interpreter in the ``.venv`` which ``uv`` creates in your
  Cleo directory. If you use a different one, replace this path with the path to your Python
  interpreter. (*hint*: if you used ``uv`` to install python for Cleo, you can find the interpreter
  path via ``uv python find``.)

* Set your Slurm options:

  replace ``<YOUR_ACCOUNT>`` and ``<YOUR_EMAIL>`` in the ``#SBATCH`` lines at the top of the job
  script with your project account and email address. You may also want to change e.g. the
  ``--time`` or ``--partition``.

Instead of editing the job script, you can also set these variables in your terminal, or add them to
your ``.bashrc`` or ``.bash_profile`` file, before using the job script, e.g.

.. code-block:: console

  export CLEO_YACYAXTROOT=your/path/to/yacyaxtroot
  export CLEO_PATH2BUILD=your/path/to/builds

The job script stops with an error if any of these paths are still set to their ``<...>``
placeholders. *Note*: ``CLEO_PATH2CLEO`` is the directory you submit the job script from, so
always submit it from your Cleo directory.

You can optionally configure the job script in the following ways:

* Choose which examples to run:

  edit the ``experiments`` list in the ``configuration`` section of the job script. Each entry
  states an example, its build configuration and its compiler, e.g. ``"as2017 serial gcc"``. You
  can instead choose one example when you submit the job script (see below).

* Choose your build configuration:

  choose which parallelism to utilise via the ``buildtype``. The options are
  ``serial``, ``threads``, ``openmp`` or ``cuda``. The default is ``openmp``. *Note*: setting ``buildtype`` to
  ``cuda`` requires you to use the GPU job script, ``gpu.sh``, and the ``gcc`` compiler.

* Choose your compiler:

  choose which compilers to use via the ``compilername``. The options are ``gcc`` or
  ``intel`` (both via MPI wrappers). *Note*: the ``cuda`` build configuration and the bubble3d
  example require you use the ``gcc`` compiler.

* Build from scratch:

  set ``CLEO_MAKE_CLEAN=true`` to delete each example's build directory before configuring Cleo
  with CMake again, e.g. after changing your compiler or build configuration.


.. _executebash_levante:

Execute the Bash Scripts
------------------------

From your Cleo directory, submit the job script to Slurm:

.. code-block:: console

  $ sbatch scripts_2/levante/cpu.sh [mode] [experiment] [buildtype] [compilername]

(or ``gpu.sh`` in place of ``cpu.sh`` to use GPUs). You must submit the job script from your
Cleo directory, because this is how the job script finds Cleo.

All the arguments are optional:

* ``mode``: ``all`` (the default) builds, compiles, runs and plots, ``build`` only builds and
  compiles, and ``run`` recompiles, runs and plots using an existing build (see
  :ref:`how the bash scripts work<bashscripts>`).

* ``experiment``: the example to run. If it is not given, every example in the job script's
  ``experiments`` list is run.

* ``buildtype`` and ``compilername``: the build configuration and compiler for the example. If
  they are not given, the defaults for Levante are used.

  *Note*: this is also true for the GPU job script, so when you choose an example with ``gpu.sh``
  also give the ``cuda`` build configuration, e.g. ``sbatch scripts_2/levante/gpu.sh all as2017 cuda``.

For example, to build Cleo and compile the executable for the Arabas and Shima 2017 example on
a login node, and then submit a job to run and plot it:

.. code-block:: console

  $ scripts_2/levante/cpu.sh build as2017
  $ sbatch scripts_2/levante/cpu.sh run as2017


The Examples
------------

.. dropdown:: Adiabatic Parcel
  :animate: fade-in

  The examples, ``as2017.py`` and ``cuspbifurc.py``, in ``examples/adiabaticparcel/`` are for a
  0-D model of a parcel of air expanding and contracting adiabatically with a two-way coupling between
  the SDM microphysics and the thermodynamics. The setup mimics that in Arabas and Shima 2017
  section 7 :cite:`arabasshima2017`. *Note*: due to numerical differences, the conditions for cusp
  bifurcation and the plots will not be exactly identical to this reference.

  .. dropdown:: a) Arabas and Shima 2017
    :animate: fade-in-slide-down

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all as2017

    The plot produced, by default called ``${CLEO_PATH2BUILD}/build_adia0d/as2017/bin/as2017fig.png``, should be
    similar to figure 5 from Arabas and Shima 2017 :cite:`arabasshima2017`.

  .. dropdown:: b) Cusp Bifurcation
    :animate: fade-in-slide-down

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all cuspbifurc

    The plots produced, by default called ``${CLEO_PATH2BUILD}/build_adia0d/cuspbifurc/bin/cuspbifurc_validation.png`` and
    ``${CLEO_PATH2BUILD}/build_adia0d/cuspbifurc/bin/cuspbifurc_SDgrowth.png`` illustrate an example of cusp bifurcation, analagous
    to the third column of figure 5 from Arabas and Shima 2017 :cite:`arabasshima2017`.


.. dropdown:: Box Model Collisions
  :animate: fade-in

  These examples, ``shima2009.py`` and ``breakup.py``, in ``examples/boxmodelcollisions/`` are for a
  0-D box model with various collision kernels. The setup mimics that in Shima et al. 2009
  section 5.1.4 :cite:`shima2009`. *Note*: due to the randomness of the initial super-droplet
  conditions and the collision algorithm, each run of these examples will not be completely identical,
  but they should be reasonably similar, and have the same mean behaviour.

  .. container:: large-text

    **The Collision Kernels:**

  *Golovin*

  The ``shima2009.py`` example models collision-coalescence using Golovin's kernel.

  The plot produced, by default called ``${CLEO_PATH2BUILD}/build_colls0d/shima2009/bin/golovin_validation.png``,
  should be similar to Fig.2(a) of Shima et al. 2009 :cite:p:`shima2009`.

  *Long*

  The ``shima2009.py`` example models collision-coalescence using Long's collision efficiency as
  given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`.

  The plot produced, by default called ``${CLEO_PATH2BUILD}/build_colls0d/shima2009/bin/long_validation_[X].png``,
  should be similar to Fig.2(b) of Shima et al. 2009 :cite:p:`shima2009`.

  *Low and List*

  The ``breakup.py`` example models collision-coalescence-rebound-breakup using the hydrodynamic
  kernel with Long's collision efficiency as given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`,
  and the coalescence/breakup/rebound probability from Low and List 1982(a) :cite:`lowlist1982a`
  (see also McFarquhar 2004 :cite:`mcfarquhar2004`). If breakup occurs, a constant
  number of fragments is produced.

  This example produces a plot, by default called ``${CLEO_PATH2BUILD}/build_colls0d/breakup/bin/lowlist_validation.png``.

  *Szakáll and Urbich*

  The ``breakup.py`` example models collision-coalescence-rebound-breakup using the hydrodynamic kernel with Long's
  collision efficiency as given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`, and the
  coalescence/breakup/rebound probability from Szakáll and Urbich 2018 :cite:`szakall2018`.
  If breakup occurs, a constant number of fragments is produced.

  This example produces a plot, by default called ``${CLEO_PATH2BUILD}/build_colls0d/breakup/bin/szakallurbich_validation.png``.

  *Testik and Straub*

  The ``breakup.py`` example models collision-coalescence-rebound-breakup using the hydrodynamic kernel with Long's
  collision efficiency as given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`, and the
  coalescence/breakup/rebound probability based on section 4 of Testik et al. 2011 (figure 12)
  :cite:`testik2011` (first proposed in :cite:`testik2009`), as well as coalescence efficiency and number of fragements
  produced from Straub et al. 2010 and Schlottke et al. 2010 respectively (:cite:`schlottke2010`, :cite:`straub2010`).

  This example produces a plot, by default called ``${CLEO_PATH2BUILD}/build_colls0d/breakup/bin/testikstraub_validation.png``.

  .. container:: large-text

    **Running the Box Model Collisions Examples:**

  .. dropdown:: a) Shima et al. 2009
    :animate: fade-in-slide-down

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all shima2009

    By default the golovin exectuable and two examples using the long executable will be compiled and
    run. You can change this by editing ``--kernels golovin long1 long2`` in the ``shima2009`` entry
    of ``scripts_2/common/experiments.sh``.

    **Golovin**

    This example models collision-coalescence using Golovin's kernel.

    The plot produced, by default called ``${CLEO_PATH2BUILD}/build_colls0d/shima2009/bin/golovin_validation.png``, should be
    comparable to Fig.2(a) of Shima et al. 2009 :cite:p:`shima2009`.

    **Long1 and Long2**

    These examples model collision-coalescence using Long's collision efficiency as given by equation
    13 of Simmel et al. 2002 :cite:`simmel2002`. The two examples use almost identical initial
    conditions and collision timesteps, as in Shima et al. 2009 :cite:p:`shima2009`.

    The plots produced, by default called ``${CLEO_PATH2BUILD}/build_colls0d/shima2009/bin/long_validation_1.png`` and
    ``${CLEO_PATH2BUILD}/build_colls0d/shima2009/bin/long_validation_2.png``, should be comparable to
    Fig.2(b) and Fig.2(c) of Shima et al. 2009 :cite:p:`shima2009`.

  .. dropdown:: b) Breakup
    :animate: fade-in-slide-down

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all breakup

    By default kernels including collision-coalescence, breakup and rebound will be compiled and
    run. You can change this by editing ``--kernels long lowlist szakallurbich testikstraub`` in the
    ``breakup`` entry of ``scripts_2/common/experiments.sh``.


.. dropdown:: Divergence Free Motion
  :animate: fade-in

  This example is run from the ``examples/divfreemotion/divfree2d.py`` script.

  1. :ref:`Configure the bash scripts<configurebash_levante>`.

  2. Submit the job script, e.g. from your Cleo directory:

  .. code-block:: console

    $ sbatch scripts_2/levante/cpu.sh all divfree2d

  This example plots the motion of super-droplets without a terminal velocity in a 2-D divergence
  free wind field. It produces a plot showing the motion of a sample of super-droplets, by default
  called ``${CLEO_PATH2BUILD}/build_divfree2d/bin/divfree2d_motion2d_validation.png``. The number of super-droplets in the domain
  should remain constant over time, as shown in the plot produced and by default called
  ``${CLEO_PATH2BUILD}/build_divfree2d/bin/divfree2d_maxnsupers_validation.png``.


.. dropdown:: 1-D Rainshafts
  :animate: fade-in

  .. dropdown:: The Original 1-D Rainshaft
    :animate: fade-in-slide-down

    This example is run from the ``examples/rainshaft1d/rainshaft1d.py`` script.

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all rainshaft1d

    Several plots and animations are produced by this example. If you would like to compare to our
    reference solutions please :ref:`contact us <contact>`.


  .. dropdown:: The EUREC4A 1-D Rainshaft
    :animate: fade-in-slide-down

    This example is a variant on the 1-d rainshaft, it runs analagously but with different inputs,
    outputs, microphysics and boundary conditions, and it produces some different plots.
    It is run from the ``examples/eurec4a1d/eurec4a1d.py`` script.

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all eurec4a1d


.. dropdown:: Constant 2-D Thermodynamics
  :animate: fade-in

  This example is run from the ``examples/constthermo2d/constthermo2d.py`` script.

  1. :ref:`Configure the bash scripts<configurebash_levante>`.

  2. Submit the job script, e.g. from your Cleo directory:

  .. code-block:: console

    $ sbatch scripts_2/levante/cpu.sh all constthermo2d

  Several plots and animations are produced by this example. If you would like to compare to our
  reference solutions please :ref:`contact us <contact>`.


.. dropdown:: 3-D Thermodynamics From File
  :animate: fade-in

  The examples, ``fromfile.py`` and ``fromfile_irreg.py``, in ``examples/fromfile/`` and
  ``examples/fromfile_irreg/`` are for a 3-D domain with time varying thermodynamics read from
  binary files. The ``fromfile_irreg.py`` example uses an irregular 3-D grid. These examples run
  the executable with MPI via ``srun``, by default with 4 MPI processes (``--ntasks=4`` in their
  entries of ``scripts_2/common/experiments.sh``).

  .. dropdown:: a) Regular Grid
    :animate: fade-in-slide-down

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all fromfile

    The plots produced, by default called ``${CLEO_PATH2BUILD}/build_fromfile/bin/ntasks4/fromfile_motion2d_validation.png``
    and ``${CLEO_PATH2BUILD}/build_fromfile/bin/ntasks4/fromfile_maxnsupers_validation.png``, show the motion of a
    sample of super-droplets and the number of super-droplets in the domain over time.

  .. dropdown:: b) Irregular Grid
    :animate: fade-in-slide-down

    1. :ref:`Configure the bash scripts<configurebash_levante>`.

    2. Submit the job script, e.g. from your Cleo directory:

    .. code-block:: console

      $ sbatch scripts_2/levante/cpu.sh all fromfile_irreg

    The plots produced are analagous to the regular grid example's, by default in
    ``${CLEO_PATH2BUILD}/build_fromfile_irreg/bin/ntasks4/``.


.. dropdown:: 3-D Bubble (Coupled to ICON via YAC)
  :animate: fade-in

  This example is run from the ``examples/bubble3d/bubble3d.py`` script. It couples Cleo to the
  dynamics of ICON's bubble test case via YAC, so it requires Cleo to be built with the ``gcc``
  compiler and :ref:`YAC and YAXT<install_yac>` to be installed. It also reads the ICON grid and data
  files given by ``orginal_icon_grid_file`` and ``orginal_icon_data_file`` in ``bubble3d.py``,
  which you may need to change to the location of these files on Levante.

  1. :ref:`Configure the bash scripts<configurebash_levante>`.

  2. Submit the job script, e.g. from your Cleo directory:

  .. code-block:: console

    $ sbatch scripts_2/levante/cpu.sh all bubble3d

  Plots of the super-droplets' motion and of the thermodynamics are produced, by default called
  ``${CLEO_PATH2BUILD}/build_bubble3d/bin/bubble_motion.png`` and ``${CLEO_PATH2BUILD}/build_bubble3d/bin/bubble_[variable].png``.


.. dropdown:: Python Bindings
  :animate: fade-in

  This example is run from the ``examples/python_bindings/python_bindings.py`` script.

  1. :ref:`Configure the bash scripts<configurebash_levante>`.

  2. Submit the job script, e.g. from your Cleo directory:

  .. code-block:: console

    $ sbatch scripts_2/levante/cpu.sh all python_bindings

  *Note*: you may have issues with python versions >= 3.14, please
  see :ref:`this note<pybind11>` for details.

  No plots are produced by this example but it should run sucessfully multiple times and produce
  ``no plotting script for python bindings example`` messages. Please note the output during
  time-stepping may not be ordered due to parallel execution.


Extension
---------
Explore ``examples/exampleplotting`` which gives examples of how to plot output from Cleo
with ``cleopy`` and ``plotcleo``, a few examples are demonstrated in the
``examples/exampleplotting/exampleplotting.py`` script.
