# AstroCoords

[![Build Status](https://github.com/HAMMERHEAD-Space/AstroCoords.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HAMMERHEAD-Space/AstroCoords.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/HAMMERHEAD-Space/AstroCoords.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/HAMMERHEAD-Space/AstroCoords.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DOI](https://zenodo.org/badge/812141681.svg)](https://doi.org/10.5281/zenodo.16954329)

AstroCoords.jl
================================

This package is intended to be a one stop shop for all things related to astrodynamics coordinate systems. In addition to being non-allocating and highly performant all transformations found here are also differentiable with compatibility with a number of different automatic and finite differencing schemas.

Currently this package implements:
- [x] Cartesian
- [x] Keplerian
- [x] Delaunay
    - [ ] Modified Delaunay 
- [x] Modified Equinoctial
- [x] Spherical
- [x] Cylindrical
- [x] Unified State Model
    - [x] USM7
    - [x] USM6
    - [x] USMEM
- [x] Milankovich
- [x] J2 Modified Equinoctial
- [ ] Generalized Modified Equinoctial Orbital Elements
- [x] EDROMO
- [x] Kustaanheimo-Stiefel
- [x] Stiefel-Scheifel

This package may eventually support Attitude Coordinates as well.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroCoords")
```

## Documentation

For more information, see the [documentation][docs-dev-url].

## Citing

If you use `AstroCoords.jl` in your work, please consider citing it.

```bibtex
@software{jordan_murphy_2025_16954330,
  author       = {Jordan Murphy},
  title        = {HAMMERHEAD-Space/AstroCoords.jl: v0.3.1},
  month        = aug,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v0.3.1},
  doi          = {10.5281/zenodo.16954330},
  url          = {https://doi.org/10.5281/zenodo.16954330},
}
```

[docs-dev-url]: https://hammerhead-space.github.io/AstroCoords.jl/dev/
[docs-stable-url]: https://hammerhead-space.github.io/AstroCoords.jl/stable/
