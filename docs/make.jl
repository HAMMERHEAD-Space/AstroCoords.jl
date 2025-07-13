using Documenter
using AstroCoords

makedocs(;
    modules=[AstroCoords],
    format=Documenter.HTML(;
        prettyurls=(!("local" in ARGS)), highlights=["yaml"], ansicolor=true
    ),
    sitename="AstroCoords.jl",
    authors="Jordan Murphy",
    pages=[
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "Coordinates" => Any[
            "coordinates/cartesian.md",
            "coordinates/delaunay.md",
            "coordinates/edromo.md",
            "coordinates/keplerian.md",
            "coordinates/milankovich.md",
            "coordinates/modEq.md",
            "coordinates/spherical.md",
            "coordinates/usm.md",
        ],
        "Utilities" => Any["util/anomalies.md", "util/quantities.md"],
        "Library" => "lib/library.md",
    ],
)

deploydocs(; repo="github.com/jmurphy6895/AstroCoords.jl.git", target="build")
