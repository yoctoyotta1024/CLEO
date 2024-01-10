.. _examples:

Examples
========

There are several test cases for CLEO, each with different builds,
domains, microphysics, coupling, Super-Droplet motion etc. They can be 
found in the ``CLEO/examples`` directory.  If you would like to
a copy of the reference solutions please :ref:`contact us <contact>`. 

Each example can be run by building CLEO with the relevant main program
and then executing the associated Python script. For running on
DKRZ's Levante HPC, there are bash scripts to help with this. 
The following instructions are intended to guide you through
running each example on Levante using their bash scripts.


.. _configurebash:

Configure the Bash Script
-------------------------

The bash script for every example follows the same layout and to use
one you will need to configure it in the following ways:

  * Use your Conda environment:

    replace the Conda environment written in the line
    stating ``source activate[…]`` with your Conda environment.

  * Use your Python version:

    replace the Python written in the line stating
    ``python=[…]`` with your Python.

  * (Optional) choose your build configuration:

    choose which thread parallelism to utilise by modifying the 
    ``kokkoshost`` and ``kokkosdevice`` flags. Please note that 
    to use CUDA parallelism you need to build and execute CLEO
    on Levante's gpu partition.

  * If you did not install CLEO in your home directory:

    modify the lines which state the ``path2CLEO`` and
    ``path2build`` to reflect this.

1) Adiabatic Parcel
-------------------
These examples are for a 0-D parcel of air expanding and
contracting adiabatically with a two-way coupling between
SDM and the thermodynamics. The setup mimics that in
Arabas and Shima 2017 section 7 :cite:`arabasshima2017`.
Note that due to numerical differences, the conditions
for cusp bifurcation and the plots will not be exactly
identical to this reference.

1. Navigate to the ``adiabaticparcel/`` directory,
e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/adiabaticparcel/


a) Arabas and Shima 2017
########################
2. :ref:`Configure<configurebash>` and execute the bash script ``as2017.sh``. 

.. code-block:: console

  $ vim as2017.sh
  $ ./as2017.sh

The plots produced, 
``~/CLEO/build/bin/as2017_fig_[x].png``, should be 
similar to the columns of figure 5 from Arabas and
Shima 2017 :cite:`arabasshima2017`.

b) Cusp Bifurcation
###################
2. :ref:`Configure<configurebash>` and execute the bash script ``cuspbifurc.sh``. 

.. code-block:: console

  $ vim cuspbifurc.sh
  $ ./cuspbifurc.sh

The plots produced, 
``~/CLEO/build/bin/cuspbifurc_validation.png`` and
``~/CLEO/build/bin/cuspbifurc_SDgrowth.png`` 
illustrate an example of cusp bifurcation, analagous to the 
third column of figure 5 from Arabas and
Shima 2017 :cite:`arabasshima2017`.

2) Box Model Collisions
-----------------------

Due to the randomness of the initial Super-Droplet conditions and
the collision-coalescence algorithm, each run of these examples
will not be completely identical, but they should be reasonably
similar, and have the same mean behaviour.

1. Navigate to the ``boxmodelcollisions/`` directory,
e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/boxmodelcollisions/

2. Configure the bash script ``shima2009.sh`` for your environment.

.. code-block:: console

  $ vim shima2009.sh

3. Execute the bash script ``shima2009.sh``. 

.. code-block:: console

  $ ./shima2009.sh

By default the golovin, long, and lowlist examples will compile
and run. You can change this by editing the arguments given to
``shima2009.py`` in the final line of the bash script.

a) Golovin
##########
This example is a 0-D box model with only collision-coalescence 
using Golovin's kernel.

The plot produced, 
``~/CLEO/build/bin/golovin_validation.png``, should be 
similar to Fig.2(a) of Shima et al. 2009 :cite:p:`shima2009`.

b) Long
#######
This example is a 0-D box model with only collision-coalescence 
using Long's collision efficiency as given by equation 13 of
Simmel et al. 2002 :cite:`simmel2002`.

The plot produced, 
``~/CLEO/build/bin/long_validation.png``, should be 
similar to Fig.2(b) of Shima et al. 2009 :cite:p:`shima2009`.

c) Low and List
###############
This example is a 0-D box model with only collision-coalescence 
using the hydrodynamic kernel with Long's collision efficiency as
given by equation 13 of Simmel et al. 2002 :cite:`simmel2002`, and the coalescence 
efficiency from Low and List 1982(a) :cite:`lowlist1982a`
(see also McFarquhar 2004 :cite:`mcfarquhar2004`).
This example produces a plot ``~/CLEO/build/bin/lowlist_validation.png``.

3) Divergence Free Motion
-------------------------

1. Navigate to the ``divfreemotion/`` directory,
e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/divfreemotion/

2. Configure the bash script ``divfree2d.sh`` for your environment.

.. code-block:: console

  $ vim divfree2d.sh

3. Execute the bash script ``divfree2d.sh``. 

.. code-block:: console

  $ ./divfree2d.sh

This example plots the motion of Super-Droplets without
sedimentation in a 2-D divergence free wind field
(see ``~/CLEO/build/bin/df2d_motion2d_validation.png``).
The number of Super-Droplets in the domain should remain
constant over time
(see ``~/CLEO/build/bin/df2d_totnsupers_validation.png``).

4) 1-D Rainshaft
------------------------------

1. Navigate to the ``rainshaft1d/`` directory,
e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/rainshaft1d/

2. Configure the bash script ``rainshaft1d.sh`` for your environment.

.. code-block:: console

  $ vim rainshaft1d.sh 

3. Execute the bash script ``rainshaft1d.sh``. 

.. code-block:: console

  $ ./rainshaft1d.sh

Several plots and animations are produced by this example. If
you would like to compare to reference solutions
please :ref:`contact us <contact>`.

5) Constant 2-D Thermodynamics 
------------------------------

1. Navigate to the ``constthermo2d/`` directory,
e.g.

.. code-block:: console

  $ cd ~/CLEO/examples/constthermo2d/

2. Configure the bash script ``constthermo2d.sh`` for your environment.

.. code-block:: console

  $ vim constthermo2d.sh 

3. Execute the bash script ``constthermo2d.sh``. 

.. code-block:: console

  $ ./constthermo2d.sh

Several plots and animations are produced by this example. If
you would like to compare to reference solutions
please :ref:`contact us <contact>`.

Extension
---------
Explore the ``exampleplotting/plotssrc`` Python module which
gives examples of how to plot output from CLEO with pySD, a few of 
which are demonstrated in the ``exampleplotting/exampleplotting.py`` 
script.