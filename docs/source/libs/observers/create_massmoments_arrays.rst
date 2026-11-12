Create Massmoments Arrays
============================

Header file: ``<libs/observers/create_massmoments_arrays.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/observers/create_massmoments_arrays.hpp>`_

Functions creating the Zarr arrays that :doc:`massmoments_observer` writes mass-moment data
into. These functions are documented alongside the rest of :doc:`massmoments_observer`'s API,
since the two headers' contents are used together and ``massmoments_observer.hpp`` already
declares them; see that page (``create_massmoment_xarray``, ``create_massmom0_xarray``,
``create_massmom1_xarray``, ``create_massmom2_xarray``) for the full signatures.

.. doxygenfile:: create_massmoments_arrays.hpp
   :project: observers
