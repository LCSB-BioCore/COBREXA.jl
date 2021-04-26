"""
    struct MATModel

Wrapper around the models loaded in dictionaries from the MATLAB representation.
"""
struct MATModel <: MetabolicModel
    mat::Dict{String,Any}
end

n_metabolites(m::MATModel)::Int = size(m.mat["S"], 1)
n_reactions(m::MATModel)::Int = size(m.mat["S"], 2)

function reactions(m::MATModel)::Vector{String}
    if haskey(m.mat, "rxns")
        reshape(m.mat["rxns"], n_reactions(m))
    else
        "rxn" .* string.(1:n_reactions(m))
    end
end

function metabolites(m::MATModel)::Vector{String}
    if haskey(m.mat, "mets")
        reshape(m.mat["mets"], n_metabolites(m))
    else
        "met" .* string.(1:n_metabolites(m))
    end
end

stoichiometry(m::MATModel) = sparse(m.mat["S"])

bounds(m::MATModel) = (
    sparse(reshape(get(m.mat, "lb", fill(-Inf, n_reactions(m), 1)), n_reactions(m))),
    sparse(reshape(get(m.mat, "ub", fill(Inf, n_reactions(m), 1)), n_reactions(m))),
)

balance(m::MATModel) =
    sparse(reshape(get(m.mat, "b", zeros(n_metabolites(m), 1)), n_metabolites(m)))

objective(m::MATModel) =
    sparse(reshape(get(m.mat, "c", zeros(n_reactions(m), 1)), n_reactions(m)))

coupling(a::MATModel)::SparseMat = spzeros(0, n_reactions(a))
n_coupling_constraints(a::MATModel)::Int = 0
coupling_bounds(a::MATModel)::Tuple{SparseVec,SparseVec} = (spzeros(0), spzeros(0))

function Base.convert(::Type{MATModel}, m::MetabolicModel)
    lb, ub = bounds(m)
    nr = n_reactions(m)
    nm = n_metabolites(m)
    return MATModel(
        Dict(
            "S" => stoichiometry(m),
            "rxns" => reshape(reactions(m), (nr, 1)),
            "mets" => reshape(metabolites(m), (nm, 1)),
            "lb" => reshape(lb, (nr, 1)),
            "ub" => reshape(ub, (nr, 1)),
            "b" => reshape(balance(m), (nm, 1)),
            "c" => reshape(objective(m), (nr, 1)),
        ),
    )
end
