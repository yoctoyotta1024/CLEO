.. _examples:

Examples
========

There are various examples of CLEO, with different build configurations, domains, microphysics,
coupling, and super-droplet motion etc. They can be found in the ``CLEO/examples`` directory. If you
would like to a copy of the reference solutions please :ref:`contact us <contact>`.

Each example can be run by building CLEO, compiling the relevant executable, and then running the
example's Python script. There are bash helper scripts for you to do all this relatively smoothly on
DKRZ's Levante HPC. The following instructions are intended to guide you through running each
example using their bash script.

Please Note: the bash script for some of the examples chooses a build configuration which uses GPUs.
To execute these scripts you will therefore need to be on a node in the GPU partition of Levante
(`see here <https://docs.dkrz.de/doc/levante/running-jobs/partitions-and-limits.html>`_
for documentation on Levante's partitions), or change the build configuration.

.. _configurebash:

Configure the Bash Scripts
--------------------------

The bash script for every example in ``scripts/levante/examples/`` provides command line
arguments to ``scripts/levante/examples/build_compile_run_plot.sh``. This script has
two steps:

1) It builds and compiles the specified exectuable(s) of CLEO by running ``scripts/levante/build_compile_cleo.sh [args]``

2) It generates input files, runs the exectuable(s), and plots the results by calling the example's Python script.


You will need to configure ``build_compile_run_plot.sh`` in the following ways:

* Use your Python version:

  replace the path in the line stating ``python=[...]`` with the path to your Python interpreter.
  (*hint*: if you used ``uv`` to install python for CLEO, you can find the interpreter path
  via ``uv python find``.)

* Set the path to your YAC and YAXT installations

  replace ``yacyaxtroot=[...]`` with the path to the directory containing your yac and yaxt
  directories, or to ``yacyaxtroot=""`` if you do not intend to run an example that requires YAC.

You can optionally configure the bash script specific to each example
(found in the same directory e.g. ``scripts/levante/examples/shima2009.sh``)
in the following ways:

* Choose your build configuration:

  choose which parallelism to utilise by modifying the ``buildtype`` parameter. The options are
  ``cuda``,  ``openmp`` or ``serial``. Note that ``buildtype="cuda"`` requires you to execute the
  script on a node in the GPU partition of Levante and may also include OpenMP parallelism.

* Choose your compiler:

  choose which compilers to use via the ``compilername`` parameter. The options are
  ``intel`` or  ``gcc`` (both via MPI wrappers). Note that the bubble3d example requires you use
  the ``gcc`` compiler.

* Choose your build directory:

  replace the path in the line stating ``path2build=[...]`` with the path you desire.

* If you did not install CLEO in your home directory:

  Ensure the lines which state the ``path2CLEO`` and ``path2build`` to reflect this.


Adiabatic Parcel
----------------
The examples, ``as2017.py`` and ``cuspbifurc.py``, in ``examples/adiabaticparcel/`` are for a
0-D model of a parcel of air expanding and contracting adiabatically with a two-way coupling between
the SDM microphysics and the thermodynamics. The setup mimics that in Arabas and Shima 2017
section 7 :cite:`arabasshima2017`. Note that due to numerical differences, the conditions for cusp
bifurcation and the plots will not be exactly identical to this reference.

a) Arabas and Shima 2017
########################

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/as2017.sh``.

2. Execute the bash script ``as2017.sh``, e.g. from your CLEO directory:

.. code-block:: console

  $ scripts/levante/examples/as2017.sh

The plot produced, by default called ``~/CLEO/build_adia0d/bin/as2017fig.png``, should be
similar to figure 5 from Arabas and Shima 2017 :cite:`arabasshima2017`.

b) Cusp Bifurcation
###################

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/cuspbifurc.sh``.

2. Execute the bash script ``cuspbifurc.sh``, e.g. from your CLEO directory:

.. code-block:: console

  $ scripts/levante/examples/cuspbifurc.sh

The plots produced, by default called ``~/CLEO/build_adia0d/bin/cuspbifurc_validation.png`` and
``~/CLEO/build_adia0d/bin/cuspbifurc_SDgrowth.png`` illustrate an example of cusp bifurcation, analagous
to the third column of figure 5 from Arabas and Shima 2017 :cite:`arabasshima2017`.


Box Model Collisions
--------------------
These examples, ``shima2009.py`` and ``breakup.py``, in ``examples/boxmodelcollisions/`` are for a
0-D box model with various collision kernels. The setup mimics that in Shima et al. 2009
section 5.1.4 :cite:`shima2009`. Note that due to the randomness of the initial super-droplet
conditions and the collision algorithm, each run of these examples will not be completely identical,
but they should be reasonably similar, and have the same mean behaviour.

The Collision Kernels
#####################

**Golovin**

The ``shima2009.py`` example models collision-coalescence using Golovin's kernel.

The plot produced, by default called ``~/CLEO/build_colls0d/[...]/bin/golovin_validation.png``,
should be similar to Fig.2(a) of Shima et al. 2009 :cite:p:`shima2009`.

**Long**

The ``shima2009.py`` example models collision-coalescence using Long's collision efficiency as
given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`.

The plot produced, by default called ``~/CLEO/build_colls0d/[...]/bin/long_validation_[X].png``,
should be similar to Fig.2(b) of Shima et al. 2009 :cite:p:`shima2009`.

**Low and List**

The ``breakup.py`` example models collision-coalescence-rebound-breakup using the hydrodynamic
kernel with Long's collision efficiency as given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`,
and the coalescence/breakup/rebound probability from Low and List 1982(a) :cite:`lowlist1982a`
(see also McFarquhar 2004 :cite:`mcfarquhar2004`). If breakup occurs, a constant
number of fragments is produced.

This example produces a plot, by default called ``~/CLEO/build_colls0d/[...]/bin/lowlist_validation.png``.

**Szakáll and Urbich**

The ``breakup.py`` example models collision-coalescence-rebound-breakup using the hydrodynamic kernel with Long's
collision efficiency as given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`, and the
coalescence/breakup/rebound probability from Szakáll and Urbich 2018 :cite:`szakall2018`.
If breakup occurs, a constant number of fragments is produced.

This example produces a plot, by default called ``~/CLEO/build_colls0d/[...]/bin/szakallurbich_validation.png``.

**Testik and Straub**

The ``breakup.py`` example models collision-coalescence-rebound-breakup using the hydrodynamic kernel with Long's
collision efficiency as given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`, and the
coalescence/breakup/rebound probability based on section 4 of Testik et al. 2011 (figure 12)
:cite:`testik2011` as well as coalescence efficiency and number of fragements produced from
Straub et al. 2010 and Schlottke et al. 2010 respectively (:cite:`schlottke2010`, :cite:`straub2010`).

This example produces a plot, by default called ``~/CLEO/build_colls0d/[...]/bin/testikstraub_validation.png``.


Running the Box Model Collisions Examples
##########################################

a) Shima et al. 2009
####################

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/shima2009.sh``.

2. Execute the bash script ``shima2009.sh``, e.g.  from your CLEO directory:

.. code-block:: console

  $ scripts/levante/examples/shima2009.sh

By default the golovin exectuable and two examples using the long executable will be compiled and
run. You can change this by editing ``script_args="[...] golovin long1 long2`` in ``shima2009.sh``.

**Golovin**

This example models collision-coalescence using Golovin's kernel.

The plot produced, by default called ``~/CLEO/build_colls0d/bin/golovin_validation.png``, should be
comparable to Fig.2(a) of Shima et al. 2009 :cite:p:`shima2009`.

**Long1 and Long2**

These examples model collision-coalescence using Long's collision efficiency as given by equation
13 of Simmel et al. 2002 :cite:`simmel2002`. The two examples use almost identical initial
conditions and collision timesteps, as in Shima et al. 2009 :cite:p:`shima2009`.

The plots produced, by default called ``~/CLEO/build_colls0d/bin/long_validation_1.png`` and
``~/CLEO/build_colls0d/bin/long_validation_2.png``, should be comparable to
Fig.2(b) and Fig.2(c) of Shima et al. 2009 :cite:p:`shima2009`.

b) Breakup
##########

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/breakup.sh``.

2. Execute the bash script ``breakup.sh``, e.g. from your CLEO directory:

.. code-block:: console

  $ scripts/levante/examples/breakup.sh

By default kernels including collision-coalescence, breakup and rebound will be compiled and
run. You can change this by editing ``script_args="[...] lowlist etc.`` in ``breakup.sh``.


Divergence Free Motion
----------------------

This example is runs from the ``examples/divfreemotion/divfree2d.py`` script.

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/divfree2d.sh``.

2. Execute the bash script ``divfree2d.sh``, e.g. from your CLEO directory:

.. code-block:: console

  $ scripts/levante/examples/divfree2d.sh

This example plots the motion of super-droplets without a terminal velocity in a 2-D divergence
free wind field. It produces a plot showing the motion of a sample of super-droplets, by default
called ``~/CLEO/build_divfree2D/bin/divfree2d_motion2d_validation.png``. The number of super-droplets in the domain
should remain constant over time, as shown in the plot produced and by default called
``~/CLEO/build_divfree2D/bin/divfree2d_totnsupers_validation.png``.


1-D Rainshaft
-------------

This example is runs from the ``examples/rainshaft1d/rainshaft1d.py`` script.

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/rainshaft1d.sh``.

2. Execute the bash script ``rainshaft1d.sh``, e.g. from your CLEO directory:

.. code-block:: console

  $ scripts/levante/examples/rainshaft1d.sh

Several plots and animations are produced by this example. If you would like to compare to our
reference solutions please :ref:`contact us <contact>`.


Constant 2-D Thermodynamics
---------------------------

This example is runs from the ``examples/constthermo2d/constthermo2d.py`` script.

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/constthermo2d.sh``

2. Execute the bash script ``constthermo2d.sh``, e.g.

.. code-block:: console

  $ scripts/levante/examples/constthermo2d.sh

Several plots and animations are produced by this example. If you would like to compare to our
reference solutions please :ref:`contact us <contact>`.


Kokkos Tools Profiling Test
---------------------------
This example, ``kokkostools.py``, in ``examples/kokkostools/`` compiles and runs the same
exectuable ``spdtest`` for four different build configurations, (1) "cuda" with CUDA and OpenMP
parallelism, (2) "openmp" with only OpenMP parallelism, (3) "threads" with only C++ threads
parallelism, and (4) "serial" without parallelism. Using the (pre-installed) Kokkos tooks'
Kernel Timer profiler, this example then outputs the time taken for each run in various ones of
CLEO's kernels.

1. :ref:`Configure the bash scripts<configurebash>`, ``scripts/levante/examples/build_compile_run_plot.sh``
and ``scripts/levante/examples/kokkostools.sh``.

2. Execute the bash script ``kokkostools.sh``, e.g.

.. code-block:: console

  $ scripts/levante/examples/kokkostools.sh

By default, a .txt file with Kokkos' simple kernel timer profiling tool data for two runs of each
of the four different build configurations is written to
``~/CLEO/build_spdtest/bin/[build_type]_[run_number]_[process_info].txt``.
The time spent in the "timestep" region can be compared with the ones
in ``~/CLEO/examples/kokkostools/spdtest_kpkerneltimer_example_solution``.

Extension
---------
Explore the ``examples/exampleplotting/plotssrc`` Python module which gives examples of how to plot output
from CLEO with ``cleopy``, a few of which are demonstrated in the ``examples/exampleplotting/exampleplotting.py``
script.
