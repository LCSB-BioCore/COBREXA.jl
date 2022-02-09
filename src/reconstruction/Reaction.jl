
"""
A small helper type for constructing reactions inline
"""
struct MetaboliteWithCoefficient
    coeff::Float64
    metabolite::Metabolite
end

const _Stoichiometry = Vector{MetaboliteWithCoefficient}
const _Stoichiometrizable = Union{Metabolite,Vector{MetaboliteWithCoefficient}}

Base.convert(::Type{_Stoichiometry}, x::Metabolite) = [MetaboliteWithCoefficient(1.0, x)]

Base.:*(a::Real, m::Metabolite) = [MetaboliteWithCoefficient(float(a), m)]

Base.:+(a::_Stoichiometrizable, b::_Stoichiometrizable) =
    vcat(convert(_Stoichiometry, a), convert(_Stoichiometry, b))

_make_reaction_dict(l, r) = Dict{String,Float64}(
    vcat(
        [mc.metabolite.id => -mc.coeff for mc in convert(_Stoichiometry, l)],
        [mc.metabolite.id => mc.coeff for mc in convert(_Stoichiometry, l)],
    ),
)

"""
    substrates → products

Make a forward-only [`Reaction`](@ref) from `substrates` and `products`.
"""
→(substrates::_Stoichiometrizable, products::_Stoichiometrizable) =
    Reaction("", _make_reaction_dict(substrates, products), :forward)

"""
    substrates ← products

Make a reverse-only [`Reaction`](@ref) from `substrates` and `products`.
"""
←(substrates::_Stoichiometrizable, products::_Stoichiometrizable) =
    Reaction("", _make_reaction_dict(substrates, products), :reverse)

"""
    substrates ↔ products

Make a bidirectional (reversible) [`Reaction`](@ref) from `substrates` and
`products`.
"""
↔(substrates::_Stoichiometrizable, products::_Stoichiometrizable) =
    Reaction("", _make_reaction_dict(substrates, products), :bidirectional)
