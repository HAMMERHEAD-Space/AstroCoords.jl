# Usage

The main way to convert between the elements is by using the AstroCoords struct. After one has been instantiated simply pass it and a gravitational parameter to a constructor of the desired element set and the package will handle the rest.

```julia
state = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

μ = 3.986004415e5

cart_state = Cartesian(state)
kep_state = Keplerian(cart_state, μ)
```

While not explicitly export if the user desired to avoid the structs, simply find the appropriate conversion inside of the coordinate_changes.jl file. Note, it make take multiple conversions to get to the desired set when using this approach.

```julia
state = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

μ = 3.986004415e5

kep_state = AstroCoords.cart2koe(state, μ)
```

## OrbitState

`OrbitState` bundles a coordinate set, an epoch, and a reference frame into a single object. The frame is encoded as a **phantom type parameter** — it appears in the type itself (e.g. `OrbitState{Cartesian{Float64}, Float64, :ICRF}`) rather than as a runtime field, enabling compile-time dispatch and zero-cost frame access.

```julia
u0 = Cartesian([-1076.2, -6765.9, -332.3, 9.357, -3.312, -1.188])
μ  = 3.986004415e5

# Construct with a plain Symbol (convenience) or Val{:ICRF} (explicit)
s = OrbitState(u0, 0.0, :ICRF)

# Accessors
coords(s)  # → the Cartesian coordinate
epoch(s)   # → 0.0
frame(s)   # → :ICRF  (reads from type, not a field)
```

### Converting coordinate type

`convert_coords` changes the coordinate set while preserving the epoch and frame:

```julia
s_kep = convert_coords(s, Keplerian, μ)
# typeof(s_kep) == OrbitState{Keplerian{Float64}, Float64, :ICRF}

# Round-trip
s_rt = convert_coords(s_kep, Cartesian, μ)
```

### Changing reference frame

`change_frame` rotates the state into a new frame. It accepts either a `FrameSystem`
(convenient, allocates a path lookup on each call) or a pre-compiled `CompiledRotation{2}`
(allocation-free, use in hot loops):

```julia
using FrameTransformations

frames = FrameSystem{2, Float64}()
add_axes!(frames, :ICRF, 1)
add_axes_rotating!(frames, :ITRF, 2, 1, my_rotation_fn)

# Convenience: Symbol + FrameSystem — allocates per call
s_itrf = change_frame(s, :ITRF, frames, μ)
# typeof(s_itrf) == OrbitState{Cartesian{Float64}, Float64, :ITRF}

# Zero-alloc: pre-compile the rotation once, reuse many times
cr = compile_rotation6(frames, :ICRF, :ITRF)
s_itrf = change_frame(s, Val{:ITRF}(), cr, μ)
```

Works with any `AstroCoord` type — the state is converted to Cartesian internally for the
rotation, then converted back to the original coordinate type.