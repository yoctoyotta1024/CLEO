Template for Observer to Write To Dataset
=========================================

Header file: ``<libs/observers/write_to_dataset_observer.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/observers/write_to_dataset_observer.hpp>`_

.. doxygenclass:: DoWriteToDataset
   :project: observers
   :private-members:
   :protected-members:
   :members:
   :undoc-members:

.. doxygenfunction:: WriteToDatasetObserver(const unsigned int interval, ParallelWriteData parallel_write)
   :project: observers

.. doxygenfunction:: WriteToDatasetObserver(const unsigned int interval, const Dataset &dataset, CollectData collect_data)
   :project: observers

.. doxygenfunction:: WriteToDatasetObserver(const unsigned int interval, const Dataset &dataset, CollectData collect_data, RaggedCount ragged_count)
   :project: observers
