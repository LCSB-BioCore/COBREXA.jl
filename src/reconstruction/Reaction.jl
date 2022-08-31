"""
$(TYPEDEF)

A small helper type for constructing reactions inline

# Fields
$(TYPEDFIELDS)
"""
struct _Stoichiometry
    s::Dict{String,Float64}
end

const _Stoichiometrizable = Union{Metabolite,_Stoichiometry}

Base.convert(::Type{_Stoichiometry}, ::Nothing) = _Stoichiometry(Dict())
Base.convert(::Type{_Stoichiometry}, m::Metabolite) = _Stoichiometry(Dict(m.id => 1.0))

Base.:*(a::Real, m::Metabolite) = _Stoichiometry(Dict(m.id => a))

"""
$(TYPEDSIGNATURES)

Shorthand for `metabolite1 + metabolite2`. Add 2 groups of [`Metabolite`](@ref)s
together to form reactions inline. Use with `+`, `*`, [`→`](@ref) and similar
operators.
"""
function Base.:+(a::_Stoichiometrizable, b::_Stoichiometrizable)
    ad = convert(_Stoichiometry, a).s
    bd = convert(_Stoichiometry, b).s
    _Stoichiometry(
        Dict(
            mid => get(ad, mid, 0.0) + get(bd, mid, 0.0) for
            mid in union(keys(ad), keys(bd))
        ),
    )
end

"""
$(TYPEDSIGNATURES)
"""
function _make_reaction_dict(r, p)
    rd = convert(_Stoichiometry, r).s
    pd = convert(_Stoichiometry, p).s
    return Dict{String,Float64}(
        mid => get(pd, mid, 0.0) - get(rd, mid, 0.0) for mid in union(keys(rd), keys(pd))
    )
end

"""
$(TYPEDSIGNATURES)

Shorthand for `substrates → products`. Make a forward-only [`Reaction`](@ref)
from `substrates` and `products`.
"""
→(substrates::Maybe{_Stoichiometrizable}, products::Maybe{_Stoichiometrizable}) =
    Reaction("", _make_reaction_dict(substrates, products), :forward)

"""
$(TYPEDSIGNATURES)

Shorthand for `substrates ← products`. Make a reverse-only [`Reaction`](@ref)
from `substrates` and `products`.
"""
←(substrates::Maybe{_Stoichiometrizable}, products::Maybe{_Stoichiometrizable}) =
    Reaction("", _make_reaction_dict(substrates, products), :reverse)

"""
$(TYPEDSIGNATURES)

Shorthand for `substrates ↔ products`. Make a bidirectional (reversible)
[`Reaction`](@ref) from `substrates` and `products`.
"""
↔(substrates::Maybe{_Stoichiometrizable}, products::Maybe{_Stoichiometrizable}) =
    Reaction("", _make_reaction_dict(substrates, products), :bidirectional)
