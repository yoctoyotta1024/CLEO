Background to the Super-Droplet Model (SDM)
===========================================

Cloud microphysics remains an integral and under-represented element of 
the climate system. Not only does this limit our understanding of clouds 
themselves, but also causes some of the largest uncertainties in climate 
modelling as a whole. Such a predicament is only exacerbated by Global 
Storm Resolving Models (GSRMs), the new generation of climate 
models which have storm-resolving resolutions O(1 km) and parametrise 
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
with Lagrangian particles, so called ‘Super-Droplets’, that act as 
representatives for the condensate populations of a cloud. 
A Super-Droplet has a multiplicity which defines how many ordinary 
condensate particles it represents. Whilst most microphysical processes, 
for example condensation and evaporation, are modelled exactly how 
ordinary condensates would be, some processes are modelled probabilistically 
instead. Indeed, the defining feature of SDM is that collisions of 
Super-Droplets are determined Monte-Carlo simulation such that the 
outcome of collisions converges towards the stochastic behaviour of a 
direct numerical simulation (DNS) as the number of Super-Droplets increases. 
As has been shown in numerous studies, SDM can reproduce the results of 
bin models at comparable computational cost, but without suffering from 
the spatial and spectral broadening caused by numerical diffusion 
:cite:`andrejczuk2010` :cite:`andrejczuk2008`
:cite:`arabasshima2013` :cite:`dziekan2019`.

In comparison with bulk and bin models, SDM has a number of important
conceptual and computational advantages :cite:`morrison2019`.