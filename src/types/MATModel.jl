"""
    struct MATModel

Wrapper around the models loaded in dictionaries from the MATLAB representation.
"""
struct MATModel <: MetabolicModel
    modeldict::Dict{String,Any}
end

n_metabolites(m::MATModel)::Int = size(m.modeldict["S"], 1)
n_reactions(m::MATModel)::Int = size(m.modeldict["S"], 2)

function reactions(m::MATModel)::Vector{String}
    if haskey(m.modeldict, "rxns")
        collect(m.modeldict["rxns"])
    else
        "rxn" .* string.(1:n_reactions(m))
    end
end

function metabolites(m::MATModel)::Vector{String}
    if haskey(m.modeldict, "mets")
        collect(m.modeldict["mets"])
    else
        "met" .* string.(1:n_metabolites(m))
    end
end

stoichiometry(m::MATModel) = sparse(m.modeldict["S"])

bounds(m::MATModel) = (
    haskey(m.modeldict, "lb") ? m.modeldict["lb"] : fill(-Inf, n_reactions(m)),
    haskey(m.modeldict, "ub") ? m.modeldict["ub"] : fill(Inf, n_reactions(m)),
)

balance(m::MATModel) =
    haskey(m.modeldict, "b") ? m.modeldict["b"] : zeros(n_metabolites(m))

objective(m::MATModel) =
    haskey(m.modeldict, "c") ? m.modeldict["c"] : zeros(n_reactions(m))

coupling(a::MATModel)::SparseMat = spzeros(0, n_reactions(a))
n_coupling_constraints(a::MATModel)::Int = 0
coupling_bounds(a::MATModel)::Tuple{SparseVec,SparseVec} = (spzeros(0), spzeros(0))
