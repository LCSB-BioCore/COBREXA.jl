"""
$(TYPEDEF)

Wrapper around the models loaded in dictionaries from the MATLAB representation.

# Fields
$(TYPEDFIELDS)
"""
struct MATModel <: AbstractMetabolicModel
    mat::Dict{String,Any}
end

Accessors.n_metabolites(m::MATModel)::Int = size(m.mat["S"], 1)
Accessors.n_reactions(m::MATModel)::Int = size(m.mat["S"], 2)

"""
$(TYPEDSIGNATURES)

Extracts reaction names from `rxns` key in the MAT file.
"""
function Accessors.reactions(m::MATModel)::Vector{String}
    if haskey(m.mat, "rxns")
        reshape(m.mat["rxns"], n_reactions(m))
    else
        "rxn" .* string.(1:n_reactions(m))
    end
end

"""
$(TYPEDSIGNATURES)

Guesses whether C in the MAT file is stored in A=[S;C].
"""
_mat_has_squashed_coupling(mat) =
    haskey(mat, "A") && haskey(mat, "b") && length(mat["b"]) == size(mat["A"], 1)


"""
$(TYPEDSIGNATURES)

Extracts metabolite names from `mets` key in the MAT file.
"""
function Accessors.metabolites(m::MATModel)::Vector{String}
    nm = n_metabolites(m)
    if haskey(m.mat, "mets")
        reshape(m.mat["mets"], length(m.mat["mets"]))[begin:nm]
    else
        "met" .* string.(1:n_metabolites(m))
    end
end

"""
$(TYPEDSIGNATURES)

Extract the stoichiometry matrix, stored under key `S`.
"""
Accessors.stoichiometry(m::MATModel) = sparse(m.mat["S"])

"""
$(TYPEDSIGNATURES)

Extracts bounds from the MAT file, saved under `lb` and `ub`.
"""
Accessors.bounds(m::MATModel) = (
    reshape(get(m.mat, "lb", fill(-Inf, n_reactions(m), 1)), n_reactions(m)),
    reshape(get(m.mat, "ub", fill(Inf, n_reactions(m), 1)), n_reactions(m)),
)

"""
$(TYPEDSIGNATURES)

Extracts balance from the MAT model, defaulting to zeroes if not present.
"""
function Accessors.balance(m::MATModel)
    b = get(m.mat, "b", spzeros(n_metabolites(m), 1))
    if _mat_has_squashed_coupling(m.mat)
        b = b[1:n_metabolites(m), :]
    end
    sparse(reshape(b, n_metabolites(m)))
end

"""
$(TYPEDSIGNATURES)

Extracts the objective from the MAT model (defaults to zeroes).
"""
Accessors.objective(m::MATModel) =
    sparse(reshape(get(m.mat, "c", zeros(n_reactions(m), 1)), n_reactions(m)))

"""
$(TYPEDSIGNATURES)

Extract coupling matrix stored, in `C` key.
"""
Accessors.coupling(m::MATModel) =
    _mat_has_squashed_coupling(m.mat) ? sparse(m.mat["A"][n_reactions(m)+1:end, :]) :
    sparse(get(m.mat, "C", zeros(0, n_reactions(m))))

"""
$(TYPEDSIGNATURES)

Extracts the coupling constraints. Currently, there are several accepted ways to store these in MATLAB models; this takes the constraints from vectors `cl` and `cu`.
"""
function Accessors.coupling_bounds(m::MATModel)
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
$(TYPEDSIGNATURES)

Extracts the possible gene list from `genes` key.
"""
function Accessors.genes(m::MATModel)
    x = get(m.mat, "genes", [])
    reshape(x, length(x))
end

"""
$(TYPEDSIGNATURES)

Extracts the associations from `grRules` key, if present.
"""
function Accessors.reaction_gene_association(m::MATModel, rid::String)
    if haskey(m.mat, "grRules")
        grr = m.mat["grRules"][findfirst(==(rid), reactions(m))]
        typeof(grr) == String ? parse_grr(grr) : nothing
    else
        nothing
    end
end

"""
$(TYPEDSIGNATURES)

Extract metabolite formula from key `metFormula` or `metFormulas`.
"""
Accessors.metabolite_formula(m::MATModel, mid::String) = maybemap(
    x -> parse_formula(x[findfirst(==(mid), metabolites(m))]),
    gets(m.mat, nothing, constants.keynames.metformulas),
)

"""
$(TYPEDSIGNATURES)

Extract metabolite charge from `metCharge` or `metCharges`.
"""
function Accessors.metabolite_charge(m::MATModel, mid::String)::Maybe{Int}
    met_charge = maybemap(
        x -> x[findfirst(==(mid), metabolites(m))],
        gets(m.mat, nothing, constants.keynames.metcharges),
    )
    maybemap(Int, isnan(met_charge) ? nothing : met_charge)
end

"""
$(TYPEDSIGNATURES)

Extract metabolite compartment from `metCompartment` or `metCompartments`.
"""
Accessors.metabolite_compartment(m::MATModel, mid::String) = maybemap(
    x -> x[findfirst(==(mid), metabolites(m))],
    gets(m.mat, nothing, constants.keynames.metcompartments),
)

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction with ID `rid`.
"""
function Accessors.reaction_stoichiometry(m::MATModel, rid::String)::Dict{String,Float64}
    ridx = first(indexin([rid], m.mat["rxns"]))
    reaction_stoichiometry(m, ridx)
end

"""
$(TYPEDSIGNATURES)

Return the stoichiometry of reaction at index `ridx`.
"""
function Accessors.reaction_stoichiometry(m::MATModel, ridx)::Dict{String,Float64}
    met_inds = findall(m.mat["S"][:, ridx] .!= 0.0)
    Dict(m.mat["mets"][met_ind] => m.mat["S"][met_ind, ridx] for met_ind in met_inds)
end

"""
$(TYPEDSIGNATURES)

Extract reaction name from `rxnNames`.
"""
Accessors.reaction_name(m::MATModel, rid::String) = maybemap(
    x -> x[findfirst(==(rid), reactions(m))],
    gets(m.mat, nothing, constants.keynames.rxnnames),
)

"""
$(TYPEDSIGNATURES)

Extract metabolite name from `metNames`.
"""
Accessors.metabolite_name(m::MATModel, mid::String) = maybemap(
    x -> x[findfirst(==(mid), metabolites(m))],
    gets(m.mat, nothing, constants.keynames.metnames),
)

# NOTE: There's no useful standard on how and where to store notes and
# annotations in MATLAB models. We therefore leave it very open for the users,
# who can easily support any annotation scheme using a custom wrapper.
# Even the (simple) assumptions about grRules, formulas and charges that we use
# here are very likely completely incompatible with >50% of the MATLAB models
# out there.

"""
$(TYPEDSIGNATURES)

Convert any metabolic model to `MATModel`.
"""
function Base.convert(::Type{MATModel}, m::AbstractMetabolicModel)
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
                default.(
                    "",
                    maybemap.(
                        x -> unparse_grr(String, x),
                        reaction_gene_association.(Ref(m), reactions(m)),
                    ),
                ),
            "metFormulas" =>
                default.(
                    "",
                    maybemap.(unparse_formula, metabolite_formula.(Ref(m), metabolites(m))),
                ),
            "metCharges" => default.(0, metabolite_charge.(Ref(m), metabolites(m))),
            "metCompartments" =>
                default.("", metabolite_compartment.(Ref(m), metabolites(m))),
        ),
    )
end
