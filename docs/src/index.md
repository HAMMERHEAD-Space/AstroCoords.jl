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

## Citing

If you use `AstroCoords.jl` in your work, please consider citing it.

```bibtex
@misc{AstroCoords.jl,
    author = {Murphy, Jordan},
    title = {AstroCoords.jl},
    year = {2024},
    publisher = {GitHub},
    journal = {GitHub repository},
    howpublished = {\\url{https://github.com/HAMMERHEAD-Space/AstroCoords.jl}}
}
```