.. _examples:

Examples
========

There are various examples of CLEO, with different build configurations, domains, microphysics,
coupling, and super-droplet motion etc. They can be found in the ``CLEO/examples`` directory. If you
would like to a copy of the reference solutions please :ref:`contact us <contact>`.

Each example can be run by building CLEO, compiling the relevant executable, and then running the
example's Python script. There are bash scripts to help you do all this on DKRZ's Levante HPC. The
following instructions are intended to guide you through running each example using their bash
script.

Please Note: the bash script for most of the examples chooses a build configuration which uses GPUs.
To execute these scripts you will therefore need to be on a node in the GPU partition of Levante
(`see here <https://docs.dkrz.de/doc/levante/running-jobs/partitions-and-limits.html>`_
for documentation on Levante's partitions), or change the build configuration.

.. _configurebash:

Configure the Bash Scripts
--------------------------

The bash script for every example provides command line arguments to ``examples/run_example.sh``. This
script has three steps:

1) It builds CLEO by running ``scripts/bash/build_cleo.sh``,

2) It compiles the specified exectuable(s) by running ``scripts/bash/compile_cleo.sh``,

3) It runs the example's Python script.


You will need to configure ```examples/run_example.sh``` in the following ways:

* Use your Conda (or Mamba) environment:

  replace the path in the line stating ``cleoenv=[…]`` with the path to your environment.

* Use your Python version:

  replace the path in the line stating ``python=[...]`` with the path to your Python interpreter.

* Set the path to your YAC and YAXT installations

  replace ``yacyaxtroot=[...]`` with the path to the directory containing your yac and yaxt
  directories, or to ``yacyaxtroot=""`` if you do not intend to run an example that requires YAC.

You can optionally configure the bash script specific to each example in the following ways:

* Choose your build configuration:

  choose which parallelism to utilise by modifying the ``buildtype`` parameter. The options are
  ``cuda``,  ``openmp`` or ``serial``. Note that ``buildtype="cuda"`` requires you to execute the script
  on a node in the GPU partition of Levante and may also include OpenMP parallelism.

* Choose your build directory:

  replace the path in the line stating ``path2build=[...]`` with the path you desire.

* If you did not install CLEO in your home directory:

  Ensure the lines which state the ``path2CLEO`` and ``path2build`` to reflect this.


Adiabatic Parcel
----------------
These examples are for a 0-D model of a parcel of air expanding and contracting adiabatically with a
two-way coupling between the SDM microphysics and the thermodynamics. The setup mimics that in
Arabas and Shima 2017 section 7 :cite:`arabasshima2017`. Note that due to numerical differences,
the conditions for cusp bifurcation and the plots will not be exactly identical to this reference.

1. Navigate to the ``adiabaticparcel/`` directory, e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/adiabaticparcel/

a) Arabas and Shima 2017
########################

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/adiabaticparcel/as2017.sh``.

3. Execute the bash script ``as2017.sh``, e.g.

.. code-block:: console

  $ ./as2017.sh

The plots produced, by default called ``~/CLEO/build_adia0D/bin/as2017fig_[x].png``, should be
similar to the columns of figure 5 from Arabas and Shima 2017 :cite:`arabasshima2017`.

b) Cusp Bifurcation
###################

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/adiabaticparcel/cuspbifurc.sh``.

3. Execute the bash script ``cuspbifurc.sh``, e.g.

.. code-block:: console

  $ ./cuspbifurc.sh

The plots produced, by default called ``~/CLEO/build_adia0D/bin/cuspbifurc_validation.png`` and
``~/CLEO/build_adia0D/bin/cuspbifurc_SDgrowth.png`` illustrate an example of cusp bifurcation, analagous
to the third column of figure 5 from Arabas and Shima 2017 :cite:`arabasshima2017`.


Box Model Collisions
--------------------
These examples are for a 0-D box model with various collision kernels. The setup mimics that in
Shima et al. 2009 section 5.1.4 :cite:`shima2009`. Note that due to the randomness of the initial
super-droplet conditions and the collision algorithm, each run of these examples will not be
completely identical, but they should be reasonably similar, and have the same mean behaviour.

The Collision Kernels
#####################

**Golovin**

This example models collision-coalescence using Golovin's kernel.

The plot produced, by default called ``~/CLEO/build_colls0D/bin/golovin_validation.png``, should be similar to
Fig.2(a) of Shima et al. 2009 :cite:p:`shima2009`.

**Long**

This example models collision-coalescence using Long's collision efficiency as given by equation
13 of Simmel et al. 2002 :cite:`simmel2002`.

The plot produced, by default called ``~/CLEO/build_colls0D/bin/long_validation_[X].png``, should be
similar to Fig.2(b) of Shima et al. 2009 :cite:p:`shima2009`.

**Low and List**

This example models collision-coalescence using the hydrodynamic kernel with Long's collision
efficiency as given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`, and the coalescence
efficiency from Low and List 1982(a) :cite:`lowlist1982a` (see also McFarquhar
2004 :cite:`mcfarquhar2004`).

This example produces a plot, by default called ``~/CLEO/build_colls0D/bin/lowlist_validation.png``.

Running the Box Model Collisions Examples
##########################################

1. Navigate to the ``boxmodelcollisions/`` directory, e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/boxmodelcollisions/

a) Shima et al. 2009
####################

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/boxmodelcollisions/shima2009.sh``.

3. Execute the bash script ``shima2009.sh``, e.g.

.. code-block:: console

  $ ./shima2009.sh

By default the golovin exectuable and two examples using the long executable will be compiled and
run. You can change this by editing ``script_args="[...] golovin long1 long2`` in ``shima2009.sh``.

**Golovin**

This example models collision-coalescence using Golovin's kernel.

The plot produced, by default called ``~/CLEO/build_colls0D/bin/golovin_validation.png``, should be
comparable to Fig.2(a) of Shima et al. 2009 :cite:p:`shima2009`.

**Long1 and Long2**

These examples model collision-coalescence using Long's collision efficiency as given by equation
13 of Simmel et al. 2002 :cite:`simmel2002`. The two examples use different initial conditions and
collision timesteps, as in Shima et al. 2009 :cite:p:`shima2009`. However the setup of the long2
example is not exactly that which makes Fig.2(c) in Shima et al. 2009.

The plots produced, by default called ``~/CLEO/build_colls0D/bin/long_validation_1.png`` and
``~/CLEO/build_colls0D/bin/long_validation_2.png``, should be comparable to
Fig.2(b) and Fig.2(c) of Shima et al. 2009 :cite:p:`shima2009`.

b) Breakup
##########

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/boxmodelcollisions/breakup.sh``.

3. Execute the bash script ``breakup.sh``, e.g.

.. code-block:: console

  $ ./breakup.sh

By default kernels including collision-coalescence, breakup and rebound will be compiled and
run. You can change this by editing ``script_args="[...] lowlist etc.`` in ``breakup.sh``.

Divergence Free Motion
----------------------

1. Navigate to the ``divfreemotion/`` directory, e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/divfreemotion/

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/boxmodelcollisions/divfree2d.sh``.

3. Execute the bash script ``divfree2d.sh``, e.g.

.. code-block:: console

  $ ./divfree2d.sh

This example plots the motion of super-droplets without a terminal velocity in a 2-D divergence
free wind field. It produces a plot showing the motion of a sample of super-droplets, by default
called ``~/CLEO/build_divfree2D/bin/df2d_motion2d_validation.png``. The number of super-droplets in the domain
should remain constant over time, as shown in the plot produced and by default called
``~/CLEO/build_divfree2D/bin/df2d_totnsupers_validation.png``.


1-D Rainshaft
-------------

1. Navigate to the ``rainshaft1d/`` directory, e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/rainshaft1d/

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/boxmodelcollisions/rainshaft1d.sh``.

3. Execute the bash script ``rainshaft1d.sh``, e.g.

.. code-block:: console

  $ ./rainshaft1d.sh

Several plots and animations are produced by this example. If you would like to compare to our
reference solutions please :ref:`contact us <contact>`.


Constant 2-D Thermodynamics
---------------------------

1. Navigate to the ``constthermo2d/`` directory, e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/constthermo2d/

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/boxmodelcollisions/constthermo2d.sh``.

3. Execute the bash script ``constthermo2d.sh``, e.g.

.. code-block:: console

  $ ./constthermo2d.sh

Several plots and animations are produced by this example. If you would like to compare to our
reference solutions please :ref:`contact us <contact>`.


Speed Test
----------
This example compiles and runs the same exectuable ``spdtest`` for three different build
configurations, (1) "cuda" with CUDA and OpenMP parallelism, (2) "openmp" with only OpenMP
parallelism, and (3) "serial" without parallelism.

1. Navigate to the ``speedtest/`` directory, e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/speedtest/

2. :ref:`Configure the bash scripts<configurebash>`, ``examples/run_example.sh`` and
``examples/boxmodelcollisions/speedtest.sh``.

3. Execute the bash script ``speedtest.sh``, e.g.

.. code-block:: console

  $ ./speedtest.sh

By default, a .txt file with the wall-clock time spent time-stepping the model for the three
different build configurations is written to ``~/CLEO/build_spdtest/bin/spd_allstats.txt``.
The times can be compared with the ones
in ``~/CLEO/examples/speedtest/speedtest_allstats_examples.txt``.

Extension
---------
Explore the ``exampleplotting/plotssrc`` Python module which gives examples of how to plot output
from CLEO with pySD, a few of which are demonstrated in the ``exampleplotting/exampleplotting.py``
script.
