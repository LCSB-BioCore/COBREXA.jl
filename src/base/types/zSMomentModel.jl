"""
    mutable struct SMomentData

Holds the already constructed SMOMENT problem.

# Fields
```
c::SparseVector{Float64, Int64}
E::SparseMatrixCSC{Float64, Int64}
d::SparseVector{Float64, Int64}
M::SparseMatrixCSC{Float64, Int64}
h::SparseVector{Float64, Int64}
reaction_map::Dict{String,Int}
metabolite_map::Dict{String,Int}
```
"""
mutable struct SMomentData
    c::SparseVector{Float64,Int64}
    E::SparseMatrixCSC{Float64,Int64}
    d::SparseVector{Float64,Int64}
    M::SparseMatrixCSC{Float64,Int64}
    h::SparseVector{Float64,Int64}
    reaction_map::Dict{String,Int}
    metabolite_map::Dict{String,Int}
end

"""
    SMomentData()

Empty constructor.
"""
SMomentData() = SMomentData(
    spzeros(0),
    spzeros(0, 0),
    spzeros(0),
    spzeros(0, 0),
    spzeros(0),
    Dict{String,Int}(),
    Dict{String,Int}(),
)

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
"""
mutable struct SMomentModel <: MetabolicModel
    smodel::StandardModel
    smomentdata::SMomentData
    enzymedata::EnzymeData
end

"""
    SMomentModel(
        model::MetabolicModel;
        reaction_kcats = Dict{String,Vector{Vector{Float64}}}(),
        reaction_protein_stoichiometry = Dict{String,Vector{Vector{Float64}}}(),
        protein_masses = Dict{String,Float64}(),
        total_protein = 0.0,
        flux_measurements = Dict{String,Tuple{Float64,Float64}}(),
    )

Construct an `SMomentModel`.

"""
function SMomentModel(
    model::MetabolicModel;
    reaction_kcats = Dict{String,Vector{Vector{Float64}}}(),
    reaction_protein_stoichiometry = Dict{String,Vector{Vector{Float64}}}(),
    protein_masses = Dict{String,Float64}(),
    total_protein = 0.0,
    flux_measurements = Dict{String,Tuple{Float64,Float64}}(),
)
    sm = convert(StandardModel, model)
    # check that input data is in correct format for smoment
    if any(length(v) > 1 for (rid, v) in reaction_kcats if has_reaction_grr(sm, rid)) ||
       any(
        length(v) > 1 for (rid, v) in reaction_protein_stoichiometry if
        haskey(reaction_kcats, rid) && has_reaction_grr(sm, rid)
    )
        @warn(
            "For SMOMENT to work correctly, no isozymes are allowed. Call `remove_slow_isozymes!` to fix the input data."
        )
    end

    smm = SMomentModel(
        sm,
        SMomentData(), # empty
        EnzymeData(
            reaction_kcats,
            reaction_protein_stoichiometry,
            protein_masses,
            total_protein;
            flux_measurements,
        ),
    )

    # build data in SMomentModel
    build_smomentmodel_internals!(smm)

    return smm
end

"""
    stoichiometry(model::SMomentModel)

Return stoichiometry matrix that includes enzymes as metabolites.
"""
function stoichiometry(model::SMomentModel)
    build_smomentmodel_internals!(model)
    return model.smomentdata.E
end

"""
    balance(model::SMomentModel)

Return stoichiometric balance.
"""
balance(model::SMomentModel) = model.smomentdata.d

"""
    objective(model::SMomentModel)

Return objective of `model`.
"""
objective(model::SMomentModel) = model.smomentdata.c

@_inherit_model_methods SMomentModel () smodel () genes
@_inherit_model_methods SMomentModel (rid::String,) smodel (rid,) reaction_gene_association reaction_stoichiometry reaction_bounds is_reaction_reversible is_reaction_forward_only is_reaction_backward_only is_reaction_unidirectional is_reaction_blocked has_reaction_isozymes has_reaction_grr

"""
    reactions(model::SMomentModel)

Returns reactions order according to stoichiometric matrix. Note, call [`genes`](@ref)
to get the order of the remaining variables.
"""
reactions(model::SMomentModel) = _order_id_to_idx_dict(model.smomentdata.reaction_map)

"""
    metabolites(model::SMomentModel)

Returns the metabolites ordered according to the stoichiometric matrix.
"""
metabolites(model::SMomentModel) = _order_id_to_idx_dict(model.smomentdata.metabolite_map)

"""
    bounds(model::SMomentModel)

Return variable bounds for `SMomentModel`.
"""
function bounds(model::SMomentModel)
    n_rxns = length(model.smomentdata.reaction_map)
    lbs = [-model.smomentdata.h[1:n_rxns]; 0]
    ubs = [model.smomentdata.h[n_rxns.+(1:n_rxns)]; model.smomentdata.h[end]]
    return lbs, ubs
end

"""
    build_smomentmodel_internals!(model::SMomentModel)

Build internal data structures used to solve SMOMENT type flux
balance analysis problems.
"""
function build_smomentmodel_internals!(model::SMomentModel)

    S, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        _build_irreversible_stoichiometric_matrix(model.smodel)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + 1

    #: equality lhs
    Se = zeros(1, n_reactions)

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(model.enzymedata.reaction_kcats, original_rid) && continue
        # these entries have kcats, only one GRR by assumption
        grr = first(reaction_gene_association(model, original_rid))
        pstoich = first(model.enzymedata.reaction_protein_stoichiometry[original_rid])
        mw = dot(pstoich, [model.enzymedata.protein_masses[gid] for gid in grr])
        kcat =
            contains(rid, "§FOR") ?
            first(model.enzymedata.reaction_kcats[original_rid])[1] :
            first(model.enzymedata.reaction_kcats[original_rid])[2]
        Se[1, col_idx] = -mw / kcat
    end

    E = [
        S zeros(n_metabolites, 1)
        Se 1.0
    ]

    #: equality rhs
    d = zeros(n_metabolites + 1)

    #: find objective
    obj_idx_orig = first(findnz(objective(model.smodel))[1])
    obj_id_orig = reactions(model.smodel)[obj_idx_orig]
    obj_id = obj_id_orig * "§FOR" # assume forward reaction is objective
    c = zeros(n_vars)
    obj_idx = reaction_map[obj_id]
    c[obj_idx] = 1.0

    #: inequality constraints
    M, h = _smoment_build_inequality_constraints(
        model,
        n_reactions,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

    #: overwrite geckomodel data
    model.smomentdata = SMomentData(
        sparse(c),
        sparse(E),
        sparse(d),
        sparse(M),
        sparse(h),
        reaction_map,
        metabolite_map,
    )

    return nothing
end

"""
    _smoment_build_inequality_constraints(
        model::SMomentModel,
        n_reactions,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

Helper function to return functions describing the inequality
constraints for smoment.
"""
function _smoment_build_inequality_constraints(
    model::SMomentModel,
    n_reactions,
    lb_fluxes,
    ub_fluxes,
    reaction_map,
)

    #: inequality lhs
    M = Array(
        [
            -I(n_reactions) zeros(n_reactions, 1)
            I(n_reactions) zeros(n_reactions, 1)
            zeros(1, n_reactions) 1
        ],
    )

    #: inequality rhs
    for original_rid in keys(model.enzymedata.flux_measurements) # only constrain if measurement available
        lb = model.enzymedata.flux_measurements[original_rid][1]
        ub = model.enzymedata.flux_measurements[original_rid][2]
        rids = [rid for rid in keys(reaction_map) if startswith(rid, original_rid)]

        if lb > 0 # forward only
            for rid in rids
                contains(rid, "§REV") && (ub_fluxes[reaction_map[rid]] = 0.0)
                contains(rid, "§FOR") &&
                    (ub_fluxes[reaction_map[rid]] = ub; lb_fluxes[reaction_map[rid]] = lb)
            end
        elseif ub < 0 # reverse only
            for rid in rids
                contains(rid, "§FOR") && (ub_fluxes[reaction_map[rid]] = 0.0)
                contains(rid, "§REV") &&
                    (ub_fluxes[reaction_map[rid]] = -lb; lb_fluxes[reaction_map[rid]] = -ub)
            end
        else # measurement does not rule our reversibility
            for rid in rids
                contains(rid, "§FOR") &&
                    (ub_fluxes[reaction_map[rid]] = ub; lb_fluxes[reaction_map[rid]] = 0)
                contains(rid, "§REV") &&
                    (ub_fluxes[reaction_map[rid]] = -lb; lb_fluxes[reaction_map[rid]] = 0)
            end
        end
    end

    h = Array([-lb_fluxes; ub_fluxes; model.enzymedata.total_protein_mass])

    return M, h
end
