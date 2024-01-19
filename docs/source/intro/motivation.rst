Motivation for CLEO
===================

In light of :ref:`such attractive properties<sdmadvatages>`, a natural
question to raise is whether SDM can be used to model warm rain more
accurately than conventional models. The discrepancies between bulk models
and observations of tropical warm rain is well documented, for example
as discussed in Schulz and Stevens 2023 :cite:`schulzstevens2023`
and vanZanten et al. 2011 :cite:`vanzanten2011`. Already SDM has been
used to improve our understanding of how precipitation formation
depends on small scale influences, for example the choice of turbulent
scheme, strength of entrainment and mixing, and the CCN size and
concentration. Such analysis uses domains with horizontal extents up to
O(10km) and resolutions no coarser than O(100m). What remains unclear is
the role that cloud microphysics has at larger scales, for example due to
the interplay between precipitation and mesoscale circulations O(100km).

By applying SDM to simulations in regional domains O(100km) with realistic
boundary conditions and large scale forcings, we would, for the first
time, have an alternative to conventional Eulerian perspectives on the
role of microphysical processes at these scales. The comparison between
such fundamentally different models gives us a new tool for assessing the
behaviour of bulk models, and showing how their use is influencing climate
simulations. Furthermore, when these simulations are weighed up with
observations, we could quantify the extent to which SDM is a better
representation not only warm rain formation, but also cloud organisation
and evolution. The recent EUREC4A campaign, has provided a multitude of
exceptionally high quality measurements which could be used for such an
assessment.

It therefore seems apparent that a new implementation of SDM is required;
capable of modelling warm rain in LES with realistic boundary
conditions and large scale forcings, and capable of application
in large regional domains, with horizontal extents O(100km). CLEO is
an attempt to build such a SDM. It strives to be a library for SDM
to model warm clouds with exceptional computational performance.
