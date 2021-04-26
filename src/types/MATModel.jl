"""
    struct MATModel

Wrapper around the models loaded in dictionaries from the MATLAB representation.
"""
struct MATModel <: MetabolicModel
    mat::Dict{String,Any}
end

n_metabolites(m::MATModel)::Int = size(m.mat["S"], 1)
n_reactions(m::MATModel)::Int = size(m.mat["S"], 2)

"""
    reactions(m::MATModel)::Vector{String}

Extracts reaction names from `rxns` key in the MAT file.
"""
function reactions(m::MATModel)::Vector{String}
    if haskey(m.mat, "rxns")
        reshape(m.mat["rxns"], n_reactions(m))
    else
        "rxn" .* string.(1:n_reactions(m))
    end
end

"""
    metabolites(m::MATModel)::Vector{String}

Extracts metabolite names from `mets` key in the MAT file.
"""
function metabolites(m::MATModel)::Vector{String}
    if haskey(m.mat, "mets")
        reshape(m.mat["mets"], n_metabolites(m))
    else
        "met" .* string.(1:n_metabolites(m))
    end
end

"""
    stoichiometry(m::MATModel)

Extract the stoichiometry matrix, stored under key `S`.
"""
stoichiometry(m::MATModel) = sparse(m.mat["S"])

"""
    bounds(m::MATModel)

Extracts bounds from the MAT file, saved under `lb` and `ub`.
"""
bounds(m::MATModel) = (
    sparse(reshape(get(m.mat, "lb", fill(-Inf, n_reactions(m), 1)), n_reactions(m))),
    sparse(reshape(get(m.mat, "ub", fill(Inf, n_reactions(m), 1)), n_reactions(m))),
)

"""
    balance(m::MATModel)

Extracts balance from the MAT model, defaulting to zeroes if not present.
"""
balance(m::MATModel) =
    sparse(reshape(get(m.mat, "b", zeros(n_metabolites(m), 1)), n_metabolites(m)))

"""
    objective(m::MATModel)

Extracts the objective from the MAT model (defaults to zeroes).
"""
objective(m::MATModel) =
    sparse(reshape(get(m.mat, "c", zeros(n_reactions(m), 1)), n_reactions(m)))

"""
    Base.convert(::Type{MATModel}, m::MetabolicModel)

Convert any metabolic model to `MATModel`.
"""
function Base.convert(::Type{MATModel}, m::MetabolicModel)
    if typeof(m) == MATModel
        return m
    end

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
