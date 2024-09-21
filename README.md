# AstroCoords

[![Build Status](https://github.com/jmurphy6895/AstroCoords.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/AstroCoords.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/jmurphy6895/AstroCoords.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/jmurphy6895/AstroForceModels.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

AstroCoords.jl
================================

This package contains various propagators of satellite trajectories for the **HAMMERHEAD.jl** ecosystem. Currently this package implements:
- [x] Cartesian
- [x] Keplerian
- [x] Delaunay
- [x] Modified Equinoctial
- [x] Spherical
- [x] Cylindrical
- [x] Unified State Model
- [x] Milankovich
- [] EDROMO
- [] Kustaanheimo-Stiefel
- [] Stiefel-Scheifel

This package may eventually support Attitude Coordinates as well.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroCoords")
```

## Documentation

For more information, see the [documentation][docs-dev-url].

[docs-dev-url]: https://jmurphy6895.github.io/AstroCoords.jl/dev/
[docs-stable-url]: https://jmurphy6895.github.io/AstroForceCoords.jl/dev/
