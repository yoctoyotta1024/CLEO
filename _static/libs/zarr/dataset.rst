Dataset
=======

Header file: ``<libs/zarr/simple_dataset.hpp>``
`[simple dataset source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/zarr/simple_dataset.hpp>`_

``SimpleDataset`` is the Zarr dataset implementation used for a run with a single MPI process
(e.g. a standalone or one-way-coupled run): every variable is written by one process, with no
collective/parallel I/O.

.. doxygenclass:: SimpleDataset
   :project: zarr
   :private-members:
   :protected-members:
   :members:
   :undoc-members:

Header file: ``<libs/zarr/collective_dataset.hpp>``
`[collective dataset source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/zarr/collective_dataset.hpp>`_

``CollectiveDataset`` is the Zarr dataset implementation used for a run split across multiple
MPI processes (e.g. a two-way-coupled run with its own domain decomposition, see
:doc:`../cartesiandomain/cartesian_decomposition`): each process writes only the data for the
gridboxes/super-droplets it owns, coordinated as a single collective write.

.. doxygenclass:: CollectiveDataset
   :project: zarr
   :private-members:
   :protected-members:
   :members:
   :undoc-members:
