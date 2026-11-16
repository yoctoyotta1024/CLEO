.. _background:

Background to Cleo
===================

This page explains the background and motivation to Cleo's creation, explained properly in the
introductions to Cleo's model description papers :cite:`bayley2026a, bayley2026b`,
and in :cite:`bayleythesis2026`.


Background to the Super-Droplet Model (SDM)
-------------------------------------------

Cloud microphysics remains an integral and under-represented element of
the climate system. Not only does this limit our understanding of clouds
themselves, but also causes some of the largest uncertainties in climate
modelling as a whole. Such a predicament is only exacerbated by Global
Storm Resolving Models (GSRMs), the new generation of climate
models which have storm-resolving resolutions O(1km) and parametrise
radiation, sub-grid turbulence, and microphysics
:cite:`slingo2022, satoh2019, stevens2019, schulthess2019`. GSRMs have
all but irradicated their parametrisations of convection, leaving
microphysical parametrisations to replace them as one of their leading
sources of uncertainty :cite:`morrison2019`. State of the art climate
models have thus accentuated the need to better understand and model
cloud microphysics.

The recently established Super-Droplet Model (SDM) is a promising
alternative to conventional bulk and bin models for cloud microphysics.
Briefly, SDM replaces traditional modelling of condensate distributions
with Lagrangian particles, so called ‘super-droplets’, that act as
representatives for the condensate populations of a cloud.
A super-droplet has a multiplicity which defines how many ordinary
condensate particles it represents. Whilst most microphysical processes,
for example condensation and evaporation, are modelled exactly how ordinary
condensates would be, collisional processes are modelled probabilistically
instead. Indeed, the defining feature of SDM is that collisions of
super-droplets are determined Monte-Carlo simulation such that the
outcome of collisions converges towards the stochastic behaviour of a
direct numerical simulation (DNS) as the number of super-droplets increases.
As has been shown in numerous studies, SDM can reproduce the results of
bin models at comparable computational cost, but without suffering from
the spatial and spectral broadening caused by numerical diffusion
:cite:`dziekan2019, arabasshima2013, andrejczuk2010, andrejczuk2008`.


.. _sdmadvatages:

The Advantages of SDM
----------------------

In comparison with bulk and bin models, SDM has a number of important
conceptual and computational advantages :cite:`morrison2019`. The Lagrangian
perspective and the reduced number of assumptions make SDM easier to
interpret physically, and therefore a more appealing model. Particularly
promising is that SDM can readily incorporate the multitude of attributes
of real droplets with only linearly increasing computational complexity.
As a consequence super-droplets can, for example, naturally convey information
about different aerosol properties, which could be used to understand how
atmospheric emissions impact rainfall or cloud reflectivity (and therefore
the Earth’s energy budget). Alternatively SDM is well suited to modelling
the plethora of cloud ice structures, as super-droplet attributes could be
designed so that different ice habits transition seamlessly between one
another. This is in stark contrast to the “curse of dimensionality” which
has plagued more complex bin and bulk models in recent decades
:cite:`grabowski2019`. The convergence properties of SDM are also a big
factor in its conceptual appeal. As the number of super-droplets increases,
SDM tends towards DNS - or in other words the closest model we have to
“true” cloud microphysics. This property arises from the mathematical
formulation of SDM and is not likewise fundamental to bulk and bin models.

From a computational perspective, whilst SDM is prohibitively expensive
for current global climate simulations, it is ideally suited to future
advances in HPC. The underlying simplicity of the model renders it highly
parallelisable and well matched to trends favouring the use of many,
extremely fast lightweight processors. As has been shown recently, this
makes SDM ideally suited to supercomputers with graphics processing units
(GPUs) :cite:`bartmanarabas2021, dziekan2019, arabas2015`.
The growth of random access memory (RAM) is also favourable for SDM.
Not only does it make SDM’s high memory usage less demanding, but it
allows for the improvement of the model’s precision by enabling more
super-sroplets to be simulated in a domain. It also allows more super-droplet
attributes and microphysical processes to be incorporated into SDM, thus
increasing the fidelity of the model to our understanding of cloud
microphysics. With regard to HPC, SDM is highly scalable and advances
in tandem.


.. _motivation:

Motivation for Cleo (a new SDM)
---------------------------------

In light of :ref:`such attractive properties<sdmadvatages>`, a natural
question to raise is whether SDM can be used to model warm rain more
accurately than conventional models. The discrepancies between bulk models
and observations of tropical warm rain is well documented, for example
as discussed in relation to the EUREC4A :cite:`schulzstevens2023`
and RICO :cite:`vanzanten2011` field campaigns. Already SDM has been
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
in large regional domains, with horizontal extents O(100km). Cleo is
an attempt to build such a SDM. It strives to be a library for SDM
to model warm clouds with exceptional computational performance.
