Background to the Super-Droplet Model (SDM)
===========================================

Cloud microphysics remains an integral and under-represented element of
the climate system. Not only does this limit our understanding of clouds
themselves, but also causes some of the largest uncertainties in climate
modelling as a whole. Such a predicament is only exacerbated by Global
Storm Resolving Models (GSRMs), the new generation of climate
models which have storm-resolving resolutions O(1km) and parametrise
radiation, sub-grid turbulence, and microphysics :cite:`slingo2022`
:cite:`satoh2019` :cite:`stevens2019` :cite:`schulthess2019`. GSRMs have
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
for example condensation and evaporation, are modelled exactly how
ordinary condensates would be, some processes are modelled probabilistically
instead. Indeed, the defining feature of SDM is that collisions of
super-droplets are determined Monte-Carlo simulation such that the
outcome of collisions converges towards the stochastic behaviour of a
direct numerical simulation (DNS) as the number of super-droplets increases.
As has been shown in numerous studies, SDM can reproduce the results of
bin models at comparable computational cost, but without suffering from
the spatial and spectral broadening caused by numerical diffusion
:cite:`dziekan2019` :cite:`arabasshima2013`
:cite:`andrejczuk2010` :cite:`andrejczuk2008`.

.. _sdmadvatages:

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
(GPUs) :cite:`bartmanarabas2021` :cite:`dziekan2019` :cite:`arabas2015`.
The growth of random access memory (RAM) is also favourable for SDM.
Not only does it make SDM’s high memory usage less demanding, but it
allows for the improvement of the model’s precision by enabling more
super-sroplets to be simulated in a domain. It also allows more super-droplet
attributes and microphysical processes to be incorporated into SDM, thus
increasing the fidelity of the model to our understanding of cloud
microphysics. With regard to HPC, SDM is highly scalable and advances
in tandem.
