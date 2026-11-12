Motion
=======

Header file: ``<libs/superdrops/motion.hpp>``
`[source] <https://github.com/yoctoyotta1024/CLEO/blob/main/libs/superdrops/motion.hpp>`_

The C++ concept for moving a super-droplet's coordinates for one (sub-)timestep, plus the
"null" (do-nothing) implementation. ``PredCorrMotion``
(:doc:`../gridboxes/predcorrmotion`) is Cleo's concrete implementation.

.. doxygenconcept:: Motion
   :project: superdrops

.. doxygenstruct:: NullMotion
   :project: superdrops
   :private-members:
   :protected-members:
   :members:
   :undoc-members:
