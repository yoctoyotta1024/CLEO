.. _examples:

Examples
========

There are several test cases for CLEO, each with different builds,
domains, microphysics, coupling, superdroplet motion etc. They can be 
found in the ``CLEO/examples`` directory. 

Each example can be run by building CLEO with the relevant main program
and then executing the associated python script. For running on
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

a) Cusp Bifurcation
###################

b) Arabas and Shima 2017
########################


2) Box Model Collisions
-----------------------

Due to the randomness of the initial superdroplet conditions and
the collision-coalescence algorithm, the plots for these examples 
will not be completely identical, but they should be reasonably
similar, especially in the mean.

1. Navigate to the ``CLEO/examples/boxmodelcollisions`` directory,
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
similar to Fig.2(a) of Shima et al. 2009.

// TODO cite properly

b) Long
#######
This example is a 0-D box model with only collision-coalescence 
using Long's collision efficiency as given by equation 13 of
Simmel et al. 2002.

The plot produced, 
``~/CLEO/build/bin/golovin_validation.png``, should be 
similar to Fig.2(b) of Shima et al. 2009.

// TODO cite properly

c) Low and List
###############
This example is a 0-D box model with only collision-coalescence 
using the hydrodynamic kernel with Long's collision efficiency as
given by equation 13 of Simmel et al. 2002, and the coalescence 
efficiency from Low and List 1982(a) (see also McFarquhar 2004).
The example produces a plot ``~/CLEO/build/bin/golovin_validation.png``.

// TODO cite properly

3) Divergence Free Motion
-------------------------


4) Constant 2-D Thermodynamics 
------------------------------