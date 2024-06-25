# Changelog

This changelog only includes the most important changes in recent updates. For a full log of all changes, please refer to git.

### Version 4.2.0
* Removed wheels. Installation is now from source to avoid conflicts with different versions of REBOUND
* Fixed gr\_full errors when Simulation was not in center of mass frame. Convergence was also improved and effect is now a bit faster.
* Improvements to unit tests (now regenerate simulationarchives each time)

### Version 4.1.0
* Fixed inneredge effect
* Updated documentation and github actions

### Version 4.0.0
* Added support for REBOUND version 4.0.0

### Version 3.12.0
* Build REBOUNDx wheels for PyPI
* Updated unit tests

### Version 3.11.0
* Added Lense-Thirring effect

### Version 3.10.0
* Added Gas dynamical friction (Generozov \& Perets 2022)

### Version 3.9.3
* Added tides\_spin for dynamical tides and spin evolution (Lu et al. 2023)
