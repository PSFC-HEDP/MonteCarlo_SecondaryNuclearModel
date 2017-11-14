MonteCarlo script that tracks particles through a user defined plasma to determine the particle's probability of undergoing a nuclear interaction. Primarily used for modeling the probability of secondary nuclear events.

MonteCarloParticleTracer.Plasma Objects are defined by radial profiles of the ion temperature, electron temperature, and mass density. They support an arbitrary number of plasma species to study the effects of impurities and they support an arbitrary number of Legendre modes to study the effects of assymetries.

Particles are sampled from MonteCarloParticleTracer.ParticleDistribution objects that define the spatial and spectral distribution of the particles. This allows the user to study the effects of different “birth distributions” (i.e. Uniform vs HotSpot models) and investigate the effects of thermal energy distributions.

This code is still very much a work in progress and has yet to be truly benchmarked against simpler models or real data.
