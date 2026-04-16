export OrbitState, convert_coords, change_frame
export coords, epoch, frame

"""
    OrbitState{C<:AstroCoord, E, F}

Bundles an orbital coordinate set with its epoch and reference frame label.
The frame `F` is encoded as a phantom type parameter (a `Symbol`), enabling
compile-time dispatch and type inference on the reference frame.

Fields:
- `coords::C` — the coordinate state (any `AstroCoord` subtype)
- `epoch::E` — epoch; accepts `Float64` (seconds past J2000 TDB) or `Tempo.Epoch`

Type parameter:
- `F` — reference frame identifier (e.g. `:ICRF`, `:ITRF`), encoded in the type
"""
struct OrbitState{C<:AstroCoord,E,F}
    coords::C
    epoch::E
end

# Primary constructor: Val{F} captures frame into type
OrbitState(coords::C, epoch::E, ::Val{F}) where {C,E,F} = OrbitState{C,E,F}(coords, epoch)

# Convenience constructor: plain Symbol — keeps call-site identical to before
function OrbitState(coords::C, epoch::E, frame::Symbol) where {C,E}
    OrbitState(coords, epoch, Val(frame))
end

# Accessors
coords(s::OrbitState) = s.coords
epoch(s::OrbitState) = s.epoch
frame(::OrbitState{C,E,F}) where {C,E,F} = F

# Epoch → seconds past J2000 conversion.
_j2000s(t::Number) = t
_j2000s(t::Epoch) = j2000s(t)

"""
    convert_coords(state::OrbitState, ::Type{T}, μ, args...) -> OrbitState

Convert the coordinate type within an `OrbitState`, preserving epoch and frame.
"""
function convert_coords(
    state::OrbitState{C,E,F}, ::Type{T}, μ::Number, args...
) where {C,E,F,T<:AstroCoord}
    return OrbitState(T(state.coords, μ, args...), state.epoch, Val(F))
end

"""
    change_frame(state::OrbitState, new_frame, frames, μ, args...) -> OrbitState

Rotate an `OrbitState` into `new_frame`. Converts to Cartesian internally, applies the
rotation, then converts back to the original coordinate type.

The `frames` argument may be a `FrameSystem` (computes the rotation at runtime via graph
lookup) or a pre-compiled `CompiledRotation{2}` (allocation-free; use `compile_rotation`
to obtain one from a `FrameSystem` when the rotation is fixed across many calls).
"""
function change_frame(
    state::OrbitState{C,E,F}, ::Val{G}, frames, μ::Number, args...
) where {C<:AstroCoord,E,F,G}
    t = _j2000s(state.epoch)
    cart = Cartesian(state.coords, μ, args...)
    R6 = rotation6(frames, F, G, t)
    cart_new = Cartesian(R6 * params(cart))
    return OrbitState(_change_frame_coords(C, cart_new, μ, args...), state.epoch, Val(G))
end

"""
    change_frame(state::OrbitState, ::Val{G}, cr::CompiledRotation{2}, μ, args...) -> OrbitState

Allocation-free `change_frame` using a pre-compiled rotation obtained via `compile_rotation`.
Bypasses the frame graph path lookup, enabling zero-allocation frame changes in hot loops.
"""
function change_frame(
    state::OrbitState{C,E,F}, ::Val{G}, cr::CompiledRotation{2}, μ::Number, args...
) where {C<:AstroCoord,E,F,G}
    t = _j2000s(state.epoch)
    cart = Cartesian(state.coords, μ, args...)
    R6 = cr(t)
    cart_new = Cartesian(R6 * params(cart))
    return OrbitState(_change_frame_coords(C, cart_new, μ, args...), state.epoch, Val(G))
end

# Cartesian → Cartesian: no coordinate conversion needed
_change_frame_coords(::Type{<:Cartesian}, cart::Cartesian, μ, args...) = cart
# Any other coord type: strip type parameters and convert from the rotated Cartesian.
# Calling C(cart, μ) on a concrete parametric type like Keplerian{Float64} would hit
# the StaticArray constructor, so we unwrap to the UnionAll first.
function _change_frame_coords(::Type{C}, cart::Cartesian, μ, args...) where {C<:AstroCoord}
    Base.typename(C).wrapper(cart, μ, args...)
end

# Convenience: accept plain Symbol for new_frame
function change_frame(state::OrbitState, new_frame::Symbol, frames, args...)
    change_frame(state, Val(new_frame), frames, args...)
end
