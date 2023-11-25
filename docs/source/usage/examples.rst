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

1) Adiabatic Parcel
-------------------
a) Cusp Bifurcation
###################

b) Arabas and Shima 2017
########################


2) Box Model Collisions
-----------------------
#. Navigate to the ``CLEO/examples/boxmodelcollisions`` directory,
e.g.
.. code-block:: console

  $ cd ~/CLEO/examples/boxmodelcollisions/

#. Configure the bash script ``shima2009.sh`` for your environment.
* Use your Conda environment:
  in the bash script ``shima2009.sh``, replace the Conda
  environment written in the line stating ``source activate[…]`` with
  your Conda environment.

* Use your Python version:
  in the bash script ``shima2009.sh``, replace the Python
  written in the line stating ``python=[…]`` with your Python.

* (Optional) choose your build configuration:
  choose which thread parallelism to utilise by modifying the 
  ``kokkoshost`` and ``kokkosdevice`` flags. Please note that 
  to use CUDA parallelism you need to build and execute CLEO
  on Levante's gpu partition.

* If you did not install CLEO in your home directory:
  change the lines stating the ``path2CLEO`` and ``path2build``
  to reflect this.

• Please note: you definitely need to make sure the python interpreter called
“python” in your conda environment is python > 3.10.4 and has the packages
matplotlib, numpy, scipy, xarray, zarr and awkward installed. They may be other
dependencies not listed here, in which case you will have to install them too.
You can install packages using the command line, e.g. via
conda activate [your conda environment]
python -m pip install [package name(s)]

a) Golovin
##########

b) Long
#######

c) Low and List
###############


3) Divergence Free Motion
-------------------------


4) Constant 2-D Thermodynamics 
------------------------------