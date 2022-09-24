
"""
$(TYPEDEF)

A "bare bones" core linear optimization problem of the form, with reaction and
metabolite names.
```
min c^T x
s.t. S x = b
      xₗ ≤ x ≤ xᵤ
```

# Fields
$(TYPEDFIELDS)
"""
mutable struct CoreModel <: MetabolicModel
    S::SparseMat
    b::SparseVec
    c::SparseVec
    xl::Vector{Float64}
    xu::Vector{Float64}
    rxns::Vector{String}
    mets::Vector{String}
    grrs::Vector{Maybe{GeneAssociation}}

    function CoreModel(
        S::MatType,
        b::VecType,
        c::VecType,
        xl::VecType,
        xu::VecType,
        rxns::StringVecType,
        mets::StringVecType,
        grrs::Vector{Maybe{GeneAssociation}} = Vector{Maybe{GeneAssociation}}(
            nothing,
            length(rxns),
        ),
    )
        all([length(b), length(mets)] .== size(S, 1)) ||
            throw(DimensionMismatch("inconsistent number of metabolites"))

        all(
            [length(c), length(xl), length(xu), length(rxns), length(grrs)] .== size(S, 2),
        ) || throw(DimensionMismatch("inconsistent number of reactions"))

        new(sparse(S), sparse(b), sparse(c), collect(xl), collect(xu), rxns, mets, grrs)
    end
end
