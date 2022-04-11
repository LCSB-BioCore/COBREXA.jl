"""
    mutable struct SMomentModel <: MetabolicModel

Construct an enzyme capacity constrained model see `Bekiaris, Pavlos Stephanos,
and Steffen Klamt. "Automatic construction of metabolic models with enzyme
constraints." BMC bioinformatics, 2020.` for implementation details.

Note, `"§"` is reserved for internal use as a delimiter, no reaction id should
contain that character. Also note, SMOMENT assumes that each reaction only has a
single enzyme (one GRR) associated with it. It is required that a model be
modified to ensure that this condition is met. For ease-of-use,
[`remove_slow_isozymes!`](@ref) is supplied to effect this. Currently only
`modifications` that change attributes of the `optimizer` are supported.

# Fields
```
reaction_ids::Vector{String}
irrev_reaction_ids::Vector{String}
metabolites::Vector{String}
c::SparseVec
S::SparseMat
b::SparseVec
xl::SparseVec
xu::SparseVec
C::SparseMat
cl::Vector{Float64}
cu::Vector{Float64}
```
"""
mutable struct SMomentModel <: MetabolicModel
    reaction_ids::Vector{String}
    irrev_reaction_ids::Vector{String}
    metabolites::Vector{String}
    c::SparseVec
    S::SparseMat
    b::SparseVec
    xl::SparseVec
    xu::SparseVec
end

"""
    stoichiometry(model::SMomentModel)

Return stoichiometry matrix that includes enzymes as metabolites.
"""
stoichiometry(model::SMomentModel) = model.S

"""
    balance(model::SMomentModel)

Return stoichiometric balance.
"""
balance(model::SMomentModel) = model.b

"""
    objective(model::SMomentModel)

Return objective of `model`.
"""
objective(model::SMomentModel) = model.c

"""
    reactions(model::SMomentModel)

Returns the reversible reactions in `model`. For 
the irreversible reactions, use [`irreversible_reactions`][@ref].
"""
reactions(model::SMomentModel) = model.reaction_ids

"""
    n_reactions(model::SMomentModel)

Returns the number of reactions in the model.
"""
n_reactions(model::SMomentModel) = length(model.reaction_ids)

"""
    irreversible_reactions(model::SMomentModel)

Returns the irreversible reactions in `model`.
"""
irreversible_reactions(model::SMomentModel) = model.irrev_reaction_ids

"""
    metabolites(model::SMomentModel)

Return the metabolites in `model`.
"""
metabolites(model::SMomentModel) = model.metabolites

"""
    n_metabolites(model::SMomentModel) = 

Return the number of metabolites in `model`.
"""
n_metabolites(model::SMomentModel) = length(metabolites(model))

"""
    bounds(model::SMomentModel)

Return variable bounds for `SMomentModel`.
"""
bounds(model::SMomentModel) = (model.xl, model.xu)

"""
    reaction_flux(model::MetabolicModel)

Helper function to get fluxes from optimization problem.
"""
function reaction_flux(model::SMomentModel)
    R = spzeros(n_reactions(model), length(model.irrev_reaction_ids) + 1)
    for (i, rid) in enumerate(reactions(model))
        for_idx = findfirst(
            x -> x == rid * "§ARM§FOR" || x == rid * "§FOR",
            model.irrev_reaction_ids,
        )
        rev_idx = findfirst(
            x -> x == rid * "§ARM§REV" || x == rid * "§REV",
            model.irrev_reaction_ids,
        )
        !isnothing(for_idx) && (R[i, for_idx] = 1.0)
        !isnothing(rev_idx) && (R[i, rev_idx] = -1.0)
    end
    return R'
end

"""
    SMomentModel(
        model::StandardModel;
        rid_isozymes = Dict{String, Vector{Isozyme}}(),
    )

Construct an `SMomentModel`.

"""
function SMomentModel(
    model::StandardModel;
    rid_isozymes = Dict{String,Vector{Isozyme}}(),
    enzyme_capacity = 0.0,
)

    # check that input data is in correct format for smoment
    if any(length(v) > 1 for v in values(rid_isozymes))
        @warn(
            "For SMOMENT to work correctly, no isozymes are allowed. Call `remove_slow_isozymes!` to fix the input data."
        )
    end

    irrevS, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        _build_irreversible_stoichiometric_matrix(model)

    #: size of resultant model
    num_reactions = size(irrevS, 2)
    num_metabolites = size(irrevS, 1)
    num_vars = num_reactions + 1

    #: equality lhs
    Se = zeros(1, num_reactions)

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(rid_isozymes, original_rid) && continue
        # these entries have kcats, only one GRR by assumption
        isozyme = first(rid_isozymes[original_rid])
        mw = sum([model.genes[gid].molar_mass * ps for (gid, ps) in isozyme.stoichiometry])
        kcat = contains(rid, "§FOR") ? first(isozyme.kcats) : last(isozyme.kcats)
        Se[1, col_idx] = -mw / kcat
    end

    S = [
        irrevS zeros(num_metabolites, 1)
        Se 1.0
    ]

    #: equality rhs
    b = zeros(num_metabolites + 1)

    #: find objective
    obj_idx_orig = first(findnz(objective(model))[1])
    obj_id_orig = reactions(model)[obj_idx_orig]
    obj_id = obj_id_orig * "§FOR" # assume forward reaction is objective
    c = spzeros(num_vars)
    obj_idx = reaction_map[obj_id]
    c[obj_idx] = 1.0

    #: bounds
    xl = sparse([lb_fluxes; 0.0])
    xu = sparse([ub_fluxes; enzyme_capacity])

    return SMomentModel(
        reactions(model),
        _order_id_to_idx_dict(reaction_map),
        _order_id_to_idx_dict(metabolite_map),
        c,
        S,
        b,
        xl,
        xu,
    )
end

"""
    change_bound!(model::SMomentModel, id; lb=nothing, ub=nothing)

Change the bound of variable in `model`. Does not change the bound if respective
bound is `nothing`. Note, for `SMomentModel`s, if the model used to construct the
`SMomentModel` has irreversible reactions, then these reactions will be
permanently irreversible in the model, i.e. changing their bounds to make them
reversible will have no effect.
"""
function change_bound!(model::SMomentModel, id; lb = nothing, ub = nothing)


    flux_for_idx =
        findfirst(x -> x == id * "§ARM§FOR" || x == id * "§FOR", model.irrev_reaction_ids)
    if !isnothing(flux_for_idx)
        if !isnothing(lb)
            if lb <= 0
                model.xl[flux_for_idx] = 0
            else
                model.xl[flux_for_idx] = lb
            end
        end
        if !isnothing(ub)
            if ub <= 0
                model.xu[flux_for_idx] = 0
            else
                model.xu[flux_for_idx] = ub
            end
        end
    end

    flux_rev_idx =
        findfirst(x -> x == id * "§ARM§REV" || x == id * "§REV", model.irrev_reaction_ids)
    if !isnothing(flux_rev_idx)
        if !isnothing(lb)
            if lb >= 0
                model.xu[flux_rev_idx] = 0
            else
                model.xu[flux_rev_idx] = -lb
            end
            if !isnothing(ub)
                if ub >= 0
                    model.xl[flux_rev_idx] = 0
                else
                    model.xl[flux_rev_idx] = -ub
                end
            end
        end
    end

    return nothing
end

"""
    change_bounds!(model::SMomentModel, ids; lbs=fill(nothing, length(ids)), ubs=fill(nothing, length(ids)))

Change the bounds of multiple variables in `model` simultaneously. See 
[`change_bound`](@ref) for details.
"""
function change_bounds!(
    model::SMomentModel,
    ids;
    lbs = fill(nothing, length(ids)),
    ubs = fill(nothing, length(ids)),
)
    for (id, lb, ub) in zip(ids, lbs, ubs)
        change_bound!(model, id; lb = lb, ub = ub)
    end
end
