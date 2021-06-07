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
    _mat_has_squashed_coupling(mat)

Guesses whether C in the MAT file is stored in A=[S;C].
"""
_mat_has_squashed_coupling(mat) =
    haskey(mat, "A") && haskey(mat, "b") && length(mat["b"]) == size(mat["A"], 1)


"""
    metabolites(m::MATModel)::Vector{String}

Extracts metabolite names from `mets` key in the MAT file.
"""
function metabolites(m::MATModel)::Vector{String}
    nm = n_metabolites(m)
    if haskey(m.mat, "mets")
        reshape(m.mat["mets"], length(m.mat["mets"]))[begin:nm]
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
function balance(m::MATModel)
    b = get(m.mat, "b", spzeros(n_metabolites(m), 1))
    if _mat_has_squashed_coupling(m.mat)
        b = b[1:n_metabolites(m), :]
    end
    sparse(reshape(b, n_metabolites(m)))
end

"""
    objective(m::MATModel)

Extracts the objective from the MAT model (defaults to zeroes).
"""
objective(m::MATModel) =
    sparse(reshape(get(m.mat, "c", zeros(n_reactions(m), 1)), n_reactions(m)))

"""
    coupling(m::MATModel)

Extract coupling matrix stored, in `C` key.
"""
coupling(m::MATModel) =
    _mat_has_squashed_coupling(m.mat) ? sparse(m.mat["A"][n_reactions(m)+1:end, :]) :
    sparse(get(m.mat, "C", zeros(0, n_reactions(m))))

"""
    coupling_bounds(m::MATModel)

Extracts the coupling constraints. Currently, there are several accepted ways to store these in MATLAB models; this takes the constraints from vectors `cl` and `cu`.
"""
function coupling_bounds(m::MATModel)
    nc = n_coupling_constraints(m)
    if _mat_has_squashed_coupling(m.mat)
        (
            sparse(fill(-Inf, nc)),
            sparse(reshape(m.mat["b"], length(m.mat["b"]))[n_reactions(m)+1:end]),
        )
    else
        (
            sparse(reshape(get(m.mat, "cl", fill(-Inf, nc, 1)), nc)),
            sparse(reshape(get(m.mat, "cu", fill(Inf, nc, 1)), nc)),
        )
    end
end

"""
    genes(m::MATModel)

Extracts the possible gene list from `genes` key.
"""
function genes(m::MATModel)
    x = get(m.mat, "genes", [])
    reshape(x, length(x))
end

"""
    reaction_gene_association(m::MATModel, rid::String)

Extracts the associations from `grRules` key, if present.
"""
function reaction_gene_association(m::MATModel, rid::String)
    if haskey(m.mat, "grRules")
        grr = m.mat["grRules"][findfirst(==(rid), reactions(m))]
        typeof(grr) == String ? _parse_grr(grr) : nothing
    else
        nothing
    end
end

"""
    metabolite_formula(m::MATModel, mid::String)

Extract metabolite formula from key `metFormula` or `metFormulas`.
"""
metabolite_formula(m::MATModel, mid::String) = _maybemap(
    x -> _parse_formula(x[findfirst(==(mid), metabolites(m))]),
    get(m.mat, "metFormula", get(m.mat, "metFormulas", nothing)),
)

"""
    metabolite_charge(m::MATModel, mid::String)

Extract metabolite charge from `metCharge` or `metCharges`.
"""
metabolite_charge(m::MATModel, mid::String) = _maybemap(
    x -> x[findfirst(==(mid), metabolites(m))],
    get(m.mat, "metCharge", get(m.mat, "metCharges", nothing)),
)

"""
    metabolite_compartment(m::MATModel, mid::String)

Extract metabolite compartment from `metCompartment` or `metCompartments`.
"""
metabolite_compartment(m::MATModel, mid::String) = _maybemap(
    x -> x[findfirst(==(mid), metabolites(m))],
    get(m.mat, "metCompartment", get(m.mat, "metCompartments", nothing)),
)


"""
    reaction_equation(model::MATModel, rxn_id::String)::Dict{String, Float64}

Return the reaction equation of reaction with id `rxn_id` in model. The reaction
equation maps metabolite ids to their stoichiometric coefficients.
"""
function reaction_equation(m::MATModel, rxn_id::String)::Dict{String, Float64}
    rxn_ind = first(indexin([rxn_id], m.mat["rxns"]))
    reaction_equation(m, rxn_ind)
end

"""
    reaction_equation(model::MATModel, rxn_ind::Int)::Dict{String, Float64}

Return the reaction equation of reaction with index `rxn_ind` in model. The reaction
equation maps metabolite ids to their stoichiometric coefficients.
"""
function reaction_equation(m::MATModel, rxn_ind::Int)::Dict{String, Float64}
    met_inds = findall(m.mat["S"][:, rxn_ind] .!= 0.0)
    Dict(m.mat["mets"][met_ind] => m.mat["S"][met_ind, rxn_ind] for met_ind in met_inds)
end

# NOTE: There's no useful standard on how and where to store notes and
# annotations in MATLAB models. We therefore leave it very open for the users,
# who can easily support any annotation scheme using a custom wrapper.
# Even the (simple) assumptions about grRules, formulas and charges that we use
# here are very likely completely incompatible with >50% of the MATLAB models
# out there.

"""
    Base.convert(::Type{MATModel}, m::MetabolicModel)

Convert any metabolic model to `MATModel`.
"""
function Base.convert(::Type{MATModel}, m::MetabolicModel)
    if typeof(m) == MATModel
        return m
    end

    lb, ub = bounds(m)
    cl, cu = coupling_bounds(m)
    nr = n_reactions(m)
    nm = n_metabolites(m)
    return MATModel(
        Dict(
            "S" => stoichiometry(m),
            "rxns" => reactions(m),
            "mets" => metabolites(m),
            "lb" => Vector(lb),
            "ub" => Vector(ub),
            "b" => Vector(balance(m)),
            "c" => Vector(objective(m)),
            "C" => coupling(m),
            "cl" => Vector(cl),
            "cu" => Vector(cu),
            "genes" => genes(m),
            "grRules" =>
                _default.(
                    "",
                    _maybemap.(
                        x -> _unparse_grr(String, x),
                        reaction_gene_association.(Ref(m), reactions(m)),
                    ),
                ),
            "metFormulas" =>
                _default.(
                    "",
                    _maybemap.(
                        _unparse_formula,
                        metabolite_formula.(Ref(m), metabolites(m)),
                    ),
                ),
            "metCharges" => _default.(0, metabolite_charge.(Ref(m), metabolites(m))),
            "metCompartments" =>
                _default.("", metabolite_compartment.(Ref(m), metabolites(m))),
        ),
    )
end
