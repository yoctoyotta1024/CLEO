Requirements
============

The following requirements ensure CLEO's build, compilation and
execution on DKRZ's Levante HPC. If they do not work,
please :ref:`contact us <contact>`.

Of course other architectures,
other compilers, versions etc. are possible, but we leave
this for you to discover.

CMake
-----
CMake minimum version 3.18.0.

On Levante it's best to use version
3.23.1 which can be loaded e.g. via the command

.. code-block:: console

  $ spack load cmake@3.23.1%gcc

Compilers
---------
A c++ compiler with the c++20 standard library is the absolute minimum.

On Levante you can use the latest gcc compilers. At the time of writing
this is gcc 11.2.0, e.g. via

.. code-block:: console

  $ module load gcc/11.2.0-gcc-11.2.0

To compile with CUDA, use Levante's latest nvhpc compilers,
e.g. via

.. code-block:: console

  $ module load nvhpc/23.7-gcc-11.2.0

Python
------
To use PySD you need Python minimum version 3.10.4.

On Levante it's a good idea to load the python3 module, e.g. via

.. code-block:: console

  $ module load python3/2022.01-gcc-11.2.0

To use PySD and to run CLEO's examples, ``matplotlib``, ``numpy``,
``scipy``, ``xarray``, ``zarr``, and ``awkward`` must be installed.
If there are other dependencies not listed here, you will have to
install them too. We kindly ask that you also
:ref:`contact us <contact>` if this is the case.

You can install packages to your conda environment e.g. via

.. code-block:: console

  $ conda activate [your conda environment]
  $ python -m pip install [package name(s)]

Pre-Commit
----------

We use pre-commit to check our code for simple issues before
submission to code review, such as before
pushing to a GitHub repository, and we reccomend you use it too.
Pre-Commit fixes for example missing semicolons, trailing whitespaces
etc., and ensures you conform to a chosen industry standard linter.

You can learn more by checking out
`pre-commit <https://pre-commit.com/>`_'s comprehensive documentation.
