CLEO Python Bindings Module
=============================

This module of Cleo contains the `pybind11 <https://pybind11.readthedocs.io>`_ C++ bindings that
expose a subset of Cleo's C++ types and functions to Python, letting Cleo be driven from, and
coupled to, Python code (e.g. a numpy-based dynamics solver, see
:doc:`coupldyn_numpy/coupldynnumpymod`). This is distinct from ``cleopy``
(:doc:`../../cleopy/cleopy`), which is pure Python with no compiled extension of its own.


Contents:
---------

.. toctree::
   :maxdepth: 1

   coupldyn_numpy/coupldynnumpymod
   cleo_python_bindings
   cleo_python_bindings_aliases
   optional_terminal_velocity
   py_cartesiandomain
   py_configuration
   py_gridboxes
   py_initialise
   py_observers
   py_runcleo
   py_superdrops
   py_zarr
