Coupldyn Numpy Module
========================

This submodule of ``cleo_python_bindings`` implements a ``CoupledDynamics`` / ``CouplingComms``
pair (see :doc:`../../runcleo/coupleddynamics` and :doc:`../../runcleo/couplingcomms`) that
exchanges gridbox state as numpy arrays, so a Python (e.g. numpy-based) dynamics solver can be
coupled to Cleo without writing any C++.


Contents:
---------

.. toctree::
   :maxdepth: 1

   numpy_dynamics
   numpy_comms
   pycoupldyn_numpy
