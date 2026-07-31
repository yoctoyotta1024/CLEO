Coalescence-Breakup-Rebound
===========================

Header file: ``<libs/superdrops/collisions/coalbure.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/superdrops/collisions/coalbure.hpp>`_

.. doxygenstruct:: DoCoalBuRe
   :project: superdrops
   :private-members:
   :protected-members:
   :members:
   :undoc-members:

.. doxygenfunction:: CoalBuRe(const unsigned int interval, const std::function<double(unsigned int)> int2realtime, const Probability collprob, const NFrags nfrags, const Flag coalbure_flag)
   :project: superdrops

.. doxygenfunction:: CoalBuRe(const unsigned int interval, const std::function<double(unsigned int)> int2realtime, const Probability collprob, const NFrags nfrags, const Flag coalbure_flag, const uint64_t seed)
   :project: superdrops
