"""
    mutable struct EnzymeData 

Holds data relevant for enzyme constrained metabolic models. 
    
Reaction turnover numbers (catalytic constants, kcats) are supplied through
`reaction_kcats`, which is a dictionary mapping reaction ids to kcats of each
isozyme. Each isozyme should have both a forward and reverse kcat, so
`reaction_kcats = Dict(rid => [[k1f, k1r], [k2f, k2r]], ...)` for `rid` with two
isozymes. The stoichiometry of each isozyme needs to be supplied by
`protein_stoichiometry`. The format is also a dictionary mapping gene ids to
their stoichiometry, e.g. `protein_stoichiometry = Dict(rid =>
[[1,1],[1,2]],...)` implies that the first isozyme of `rid` is composed of two
subunits, each present once in the protein, while the second isozyme is composed
of two subunits, but the second subunit is present twice in the isozyme. The
order of each entry in `reaction_kcats` and `reaction_protein_stoichiometry` is
taken to be the same as the order returned when calling
[`reaction_gene_association`](@ref) on the model. The protein masses (in molar
mass units) for each gene in the model should be supplied through
`protein_masses`. The format is a dictionary of gene ids mapped to molar masses. 

Total enzyme capacity (sum of all enzyme concentrations multiplied by their
molar mass) is constrained by `total_protein_mass`, a unitless mass fraction of
enzyme mass to cell dry mass. The reaction fluxes and protein concentrations can
be bounded by `flux_measurements` and `protein_measurements` respectively. Both
lower and upper bounds need to be supplied (as a tuple) if a reaction flux is to
be bounded, likewise with protein concentration bounds. 

# Fields
```    
reaction_kcats::Dict{String,Vector{Vector{Float64}}} # rid => [[for, rev], ...]
reaction_protein_stoichiometry::Dict{String,Vector{Vector{Float64}}} # rid => [[stoich, stoich,...], ...]
protein_masses::Dict{String,Float64}
total_protein_mass::Float64
flux_measurements::Dict{String,Tuple{Float64,Float64}} # rid => (lb, ub)
protein_measurements::Dict{String,Tuple{Float64,Float64}} # pid => (lb, ub)
```
"""
mutable struct EnzymeData 
    reaction_kcats::Dict{String,Vector{Vector{Float64}}} # rid => [[for, rev], ...]
    reaction_protein_stoichiometry::Dict{String,Vector{Vector{Float64}}} # rid => [[stoich, stoich,...], ...]
    protein_masses::Dict{String,Float64}
    total_protein_mass::Float64
    flux_measurements::Dict{String,Tuple{Float64,Float64}} # rid => (lb, ub)
    protein_measurements::Dict{String,Tuple{Float64,Float64}} # pid => (lb, ub)
end

"""
    EnzymeData(
        reaction_kcats,
        reaction_protein_stoichiometry,
        protein_masses,
        total_protein;
        flux_measurements = Dict{String,Tuple{Float64,Float64}}(),
        protein_measurements = Dict{String,Tuple{Float64,Float64}}(),
    )

Constructor for `EnzymeData`.
"""
EnzymeData(
    reaction_kcats,
    reaction_protein_stoichiometry,
    protein_masses,
    total_protein;
    flux_measurements = Dict{String,Tuple{Float64,Float64}}(),
    protein_measurements = Dict{String,Tuple{Float64,Float64}}(),
) = EnzymeData(
    reaction_kcats,
    reaction_protein_stoichiometry,
    protein_masses,
    total_protein,
    flux_measurements,
    protein_measurements,
)

"""
    mutable struct GeckoData

Holds the already constructed GECKO problem.

# Fields
```    
c::SparseVector{Float64, Int64}
E::SparseMatrixCSC{Float64, Int64}
d::SparseVector{Float64, Int64}
M::SparseMatrixCSC{Float64, Int64}
h::SparseVector{Float64, Int64}
reaction_map::Dict{String,Int}
metabolite_map::Dict{String,Int}
protein_ids::Vector{String}
```
"""
mutable struct GeckoData
    c::SparseVector{Float64,Int64}
    E::SparseMatrixCSC{Float64,Int64}
    d::SparseVector{Float64,Int64}
    M::SparseMatrixCSC{Float64,Int64}
    h::SparseVector{Float64,Int64}
    reaction_map::Dict{String,Int}
    metabolite_map::Dict{String,Int}
    protein_ids::Vector{String}
end

"""
    GeckoData()

Empty constructor.
"""
GeckoData() = GeckoData(
    spzeros(0),
    spzeros(0, 0),
    spzeros(0),
    spzeros(0, 0),
    spzeros(0),
    Dict{String,Int}(),
    Dict{String,Int}(),
    Vector{String}(),
)

"""
    mutable struct GeckoModel <: MetabolicModel

A model that incorporates enzyme capacity and kinetic constraints via the GECKO
formulation. See `Sánchez, Benjamín J., et al. "Improving the phenotype
predictions of a yeast genome‐scale metabolic model by incorporating enzymatic
constraints." Molecular systems biology, 2017.` for implementation details.

Note, since the model uses irreversible reactions internally, `"§FOR"` (for the
forward direction) and `"§REV"` (for the reverse direction) is appended to each
reaction internally. Futhermore, `"§"` is reserved for internal use as a
delimiter, no reaction id should contain that character. The units depend on
those used in `enzymedata.reaction_kcats` and `enzymedata.protein_masses`. Only
the protein and reaction flux bounds are optional parameters, all other
parameters must be supplied to the `enzymedata` field. Only reactions with kcats
will have enzyme bounds associated with them, but all isozymes are assumed to
have data if data is supplied. Currently only `modifications` that change
attributes of the `optimizer` are supported.

To actually run GECKO, call [`flux_balance_analysis`](@ref) on a `GeckoModel` 
to run an analysis on it.

See also: [`StandardModel`](@ref)

# Fields
```
smodel::StandardModel
geckodata::GeckoData
enzymedata::EnzymeData
```
"""
mutable struct GeckoModel <: MetabolicModel
    smodel::StandardModel
    geckodata::GeckoData
    enzymedata::EnzymeData
end

"""
    GeckoModel(
        model::MetabolicModel;
        reaction_kcats = Dict{String,Vector{Vector{Float64}}}(),
        reaction_protein_stoichiometry = Dict{String,Vector{Vector{Float64}}}(),
        protein_masses = Dict{String,Float64}(),
        total_protein = 0.0,
        flux_measurements = Dict{String,Tuple{Float64,Float64}}(),
        protein_measurements = Dict{String,Tuple{Float64,Float64}}(),
    )

Constructor for `GeckoModel`. 
"""
function GeckoModel(
    model::MetabolicModel;
    reaction_kcats = Dict{String,Vector{Vector{Float64}}}(),
    reaction_protein_stoichiometry = Dict{String,Vector{Vector{Float64}}}(),
    protein_masses = Dict{String,Float64}(),
    total_protein = 0.0,
    flux_measurements = Dict{String,Tuple{Float64,Float64}}(),
    protein_measurements = Dict{String,Tuple{Float64,Float64}}(),
)
    gm = GeckoModel(
        convert(StandardModel, model),
        GeckoData(), # empty
        EnzymeData(
            reaction_kcats,
            reaction_protein_stoichiometry,
            protein_masses,
            total_protein;
            flux_measurements,
            protein_measurements,
        ),
    )

    # build data in GeckoModel
    build_geckomodel_internals!(gm)

    return gm
end

"""
    stoichiometry(model::GeckoModel)

Return stoichiometry matrix that includes enzymes as metabolites.
"""
function stoichiometry(model::GeckoModel)
    build_geckomodel_internals!(model)
    return model.geckodata.E
end

"""
    balance(model::GeckoModel)

Return stoichiometric balance.
"""
balance(model::GeckoModel) = model.geckodata.d

"""
    objective(model::GeckoModel)

Return objective of `model`.
"""
objective(model::GeckoModel) = model.geckodata.c

@_inherit_model_methods GeckoModel (rid::String,) smodel (rid,) reaction_gene_association reaction_stoichiometry reaction_bounds is_reaction_reversible is_reaction_forward_only is_reaction_backward_only is_reaction_unidirectional is_reaction_blocked has_reaction_isozymes has_reaction_grr

"""
    reactions(model::GeckoModel)

Returns reactions order according to stoichiometric matrix. Note, call [`genes`](@ref)
to get the order of the remaining variables.
"""
reactions(model::GeckoModel) = _order_id_to_idx_dict(model.geckodata.reaction_map)

"""
    metabolites(model::GeckoModel)

Returns the metabolites ordered according to the stoichiometric matrix.
"""
metabolites(model::GeckoModel) = _order_id_to_idx_dict(model.geckodata.metabolite_map)

"""
    genes(model::GeckoModel)

Returns the genes (proteins) in the order as they appear as variables in the model.
"""
genes(model::GeckoModel) = model.geckodata.protein_ids

"""
    _order_id_to_idx_dict(id_to_idx_dict)

Return the keys of `id_to_idx_dict` sorted by the values, which
are taken to be the indices. This is a helper function for 
[`reactions`](@ref) and [`metabolites`](@ref).
"""
function _order_id_to_idx_dict(dmap)
    ks = collect(keys(dmap))
    vs = collect(values(dmap))
    return ks[sortperm(vs)]
end

"""
    bounds(model::GeckoModel)

Return variable bounds for `GeckoModel`.
"""
function bounds(model::GeckoModel)
    n_rxns = length(model.geckodata.reaction_map)
    n_prots = length(model.geckodata.protein_ids)
    lbs = [-model.geckodata.h[1:n_rxns]; -model.geckodata.h[2*n_rxns.+(1:n_prots)]]
    ubs = [model.geckodata.h[n_rxns.+(1:n_rxns)]; model.geckodata.h[2*n_rxns+n_prots.+(1:n_prots)]]
    return lbs, ubs
end

"""
    enzyme_capacity(model::GeckoModel)

Return enzyme capacity inequality constraint vector and bound, or nothing 
if it doesn't exist in the model.
"""
enzyme_capacity(model::GeckoModel) = (model.geckodata.M[end, :], model.geckodata.h[end])

"""
    build_geckomodel_internals!(model::GeckoModel)

Lower level function that updates the matrix form of a model with enzyme
capacity constraints, in GECKO format.

Specifically, updates `model.geckodata` with the vector and matrix coefficients `c,
E, d, M, h` satisfying 
```
opt cᵀ * x
s.t.    E * x = d 
        M * x ≤ h
```
as well as `reaction_map, metabolite_map, protein_ids`, where 
`reaction_map` shows the order of the columns (reactions) in `E`. Proteins 
are ordered according to `protein_ids`, and follow after reactions.
"""
function build_geckomodel_internals!(model::GeckoModel)
    S, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        _build_irreversible_stoichiometric_matrix(model.smodel)

    #: find all gene products that have kcats associated with them
    protein_ids = _get_proteins_with_kcats(model)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_proteins = length(protein_ids)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + n_proteins

    #: equality lhs
    E_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeffs = Vector{Float64}(),
    )

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        contains(rid, "§ARM") && continue
        !haskey(model.enzymedata.reaction_kcats, original_rid) && continue

        # these entries have kcats
        if contains(rid, "§ISO")
            iso_num = parse(
                Int,
                replace(
                    first(filter(startswith("ISO"), split(rid, "§")[2:end])),
                    "ISO" => "",
                ),
            )
        else # only one enzyme
            iso_num = 1
        end

        # add all entries to column of matrix
        _add_enzyme_variable(
            model,
            iso_num, # only one enzyme
            rid,
            original_rid,
            E_components,
            col_idx,
            protein_ids,
        )
    end

    Se = sparse(
        E_components.row_idxs,
        E_components.col_idxs,
        E_components.coeffs,
        n_proteins,
        n_reactions,
    )

    E = [
        S zeros(n_metabolites, n_proteins)
        Se I(n_proteins)
    ]

    #: equality rhs
    d = zeros(n_metabolites + n_proteins)

    #: find objective
    obj_idx_orig = first(findnz(objective(model.smodel))[1])
    obj_id_orig = reactions(model.smodel)[obj_idx_orig]
    obj_id = obj_id_orig * "§FOR" # assume forward reaction is objective
    c = zeros(n_vars)
    obj_idx = reaction_map[obj_id]
    c[obj_idx] = 1.0

    #: inequality constraints
    M, h = _gecko_build_inequality_constraints(
        model,
        protein_ids,
        n_reactions,
        n_proteins,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

    #: overwrite geckomodel data
    model.geckodata = GeckoData(
        sparse(c),
        sparse(E),
        sparse(d),
        sparse(M),
        sparse(h),
        reaction_map,
        metabolite_map,
        protein_ids,
    )

    return nothing
end

"""
    _gecko_build_inequality_constraints(
        model::GeckoModel,
        protein_ids,
        n_reactions,
        n_proteins,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

Helper function to build inequality constraints. Returns the inequality constraint in matrix format.
"""
function _gecko_build_inequality_constraints(
    model::GeckoModel,
    protein_ids,
    n_reactions,
    n_proteins,
    lb_fluxes,
    ub_fluxes,
    reaction_map,
)
    #: inequality lhs
    mw_proteins = [model.enzymedata.protein_masses[pid] for pid in protein_ids]
    M = Array(
        [
            -I(n_reactions) zeros(n_reactions, n_proteins)
            I(n_reactions) zeros(n_reactions, n_proteins)
            zeros(n_proteins, n_reactions) -I(n_proteins)
            zeros(n_proteins, n_reactions) I(n_proteins)
            zeros(1, n_reactions) mw_proteins'
        ],
    )

    #: inequality rhs
    for original_rid in keys(model.enzymedata.flux_measurements) # only constrain if measurement available
        lb = model.enzymedata.flux_measurements[original_rid][1]
        ub = model.enzymedata.flux_measurements[original_rid][2]
        rids = [rid for rid in keys(reaction_map) if startswith(rid, original_rid)]
        filter!(x -> !contains(x, "§ISO"), rids) # remove isozyme partial reactions (ARM reactions take care of these)

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

    lb_proteins = [
        haskey(model.enzymedata.protein_measurements, pid) ? model.enzymedata.protein_measurements[pid][1] : 0.0 for pid in protein_ids
    ]
    ub_proteins = [
        haskey(model.enzymedata.protein_measurements, pid) ? model.enzymedata.protein_measurements[pid][2] :
        1000.0 for pid in protein_ids
    ]

    h = Array([-lb_fluxes; ub_fluxes; -lb_proteins; ub_proteins; model.enzymedata.total_protein_mass])

    return M, h
end

"""
    _build_irreversible_stoichiometric_matrix(model::StandardModel)

Return the stoichiometric matrix. All reactions are forward only i.e. only 
positive fluxes are allowed. Include arm reactions.
"""
function _build_irreversible_stoichiometric_matrix(model::StandardModel)
    # components used to build stoichiometric matrix
    S_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeffs = Vector{Float64}(),
        lbs = Vector{Float64}(),
        ubs = Vector{Float64}(),
    )

    # establish the ordering in a named tuple
    idxs = ( #: pseudo metabolites and reactions are added to model
        met_idxs = Dict{String,Int}(),
        rxn_idxs = Dict{String,Int}(),
        max_rxn_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
        max_met_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
        pseudo_met_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
    )
    #TODO for the counter thing, basically I wanted e.g. max_rxn_idx = 1 and then update it, 
    #TODO but named tuples are immutable... :(

    # fill the matrix entries
    #: blocked treated as reversible because unclear what direction the reaction would go
    for rid in reactions(model)
        if has_reaction_grr(model, rid) && has_reaction_isozymes(model, rid)
            if is_reaction_unidirectional(model, rid)
                dir = is_reaction_forward_only(model, rid) ? "§FOR" : "§REV"
                _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
            elseif is_reaction_reversible(model, rid) || is_reaction_blocked(model, rid)
                _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§FOR")
                _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§REV")
            else
                @warn "Unhandled bound type for $rid"
            end
        else # no grr or single enzyme only
            if is_reaction_unidirectional(model, rid)
                dir = is_reaction_forward_only(model, rid) ? "§FOR" : "§REV"
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
            elseif is_reaction_reversible(model, rid) || is_reaction_blocked(model, rid)
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§FOR")
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§REV")
            else
                @warn "Unhandled bound type for $rid"
            end
        end
    end

    S = sparse(
        S_components.row_idxs,
        S_components.col_idxs,
        S_components.coeffs,
        length(idxs.met_idxs),
        length(idxs.rxn_idxs),
    )

    return S, S_components.lbs, S_components.ubs, idxs.rxn_idxs, idxs.met_idxs
end

"""
    _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)

Add entries to the components that will be used to build the stoichiometric
matrix. Simple variant that does not deal with isozymes and arm reactions.
"""
function _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
    idxs.rxn_idxs[rid*dir] = idxs.max_rxn_idx[1]
    idxs.max_rxn_idx[1] += 1
    fix_sign = dir == "§FOR" ? 1 : -1 # change direction of reaction
    for (mid, coeff) in reaction_stoichiometry(model, rid)
        if !haskey(idxs.met_idxs, mid)
            idxs.met_idxs[mid] = idxs.max_met_idx[1]
            idxs.max_met_idx[1] += 1
        end
        push!(S_components.row_idxs, idxs.met_idxs[mid])
        push!(S_components.col_idxs, idxs.rxn_idxs[rid*dir])
        push!(S_components.coeffs, fix_sign * coeff)
    end
    lb, ub = abs.(reaction_bounds(model, rid)) # assumes lb < ub
    if dir == "§FOR"
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, lb)
        push!(S_components.ubs, ub)
    else
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, ub)
        push!(S_components.ubs, lb)
    end
end

"""
    _add_isozyme_to_irrev_stoich_mat(model::GeckoModel, rid, idxs, S_components, dir)

Add entries to the components that will be used to build the stoichiometric matrix.
Complex variant that deals with isozymes and arm reactions.
"""
function _add_isozyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
    # add pseudo metabolite
    pm = "§PM$(idxs.pseudo_met_idx[1])"
    idxs.pseudo_met_idx[1] += 1
    idxs.met_idxs[pm] = idxs.max_met_idx[1]
    idxs.max_met_idx[1] += 1
    # find half reactions to get arm reaction
    lhs = []
    rhs = []
    for (mid, coeff) in reaction_stoichiometry(model, rid)
        if !haskey(idxs.met_idxs, mid)
            idxs.met_idxs[mid] = idxs.max_met_idx[1]
            idxs.max_met_idx[1] += 1
        end
        if coeff <= 0
            push!(lhs, (mid, coeff))
        else
            push!(rhs, (mid, coeff))
        end
    end
    product_half_reaction = dir == "§FOR" ? rhs : lhs
    reagent_half_reaction = dir == "§FOR" ? lhs : rhs
    # add arm reaction
    fix_sign = dir == "§FOR" ? 1 : -1 # change direction of reaction
    pr = rid * "§ARM" * dir
    idxs.rxn_idxs[pr] = idxs.max_rxn_idx[1] #! this needs to get added first because of blocked possibility
    idxs.max_rxn_idx[1] += 1
    push!(S_components.row_idxs, idxs.met_idxs[pm])
    push!(S_components.col_idxs, idxs.rxn_idxs[pr])
    push!(S_components.coeffs, 1)
    for (mid, coeff) in reagent_half_reaction
        push!(S_components.row_idxs, idxs.met_idxs[mid])
        push!(S_components.col_idxs, idxs.rxn_idxs[pr])
        push!(S_components.coeffs, fix_sign * coeff)
    end
    # add bounds for ARM reaction that corresponds to original model's bounds
    lb, ub = abs.(reaction_bounds(model, rid)) # assumes lb < ub
    if dir == "§FOR"
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, lb)
        push!(S_components.ubs, ub)
    else
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, ub)
        push!(S_components.ubs, lb)
    end
    # add isozyme reactions
    for (i, _) in enumerate(reaction_gene_association(model, rid))
        iso_rid = rid * "§ISO$i" * dir
        idxs.rxn_idxs[iso_rid] = idxs.max_rxn_idx[1]
        idxs.max_rxn_idx[1] += 1
        push!(S_components.row_idxs, idxs.met_idxs[pm])
        push!(S_components.col_idxs, idxs.rxn_idxs[iso_rid])
        push!(S_components.coeffs, -1)
        for (mid, coeff) in product_half_reaction
            push!(S_components.row_idxs, idxs.met_idxs[mid])
            push!(S_components.col_idxs, idxs.rxn_idxs[iso_rid])
            push!(S_components.coeffs, fix_sign * coeff)
        end
        # add bounds
        push!(S_components.lbs, 0)
        if is_reaction_blocked(model, rid)
            push!(S_components.ubs, 0)
        else
            push!(S_components.ubs, 1000) # arbitrary upper bound
        end
    end
end

"""
    _add_enzyme_variable(
        model::GeckoModel,
        iso_num,
        rid,
        original_rid,
        E_components,
        col_idx,
        protein_ids,
    )

Helper function to add an column into the enzyme stoichiometric matrix.
"""
function _add_enzyme_variable(
    model::GeckoModel,
    iso_num,
    rid,
    original_rid,
    E_components,
    col_idx,
    protein_ids,
)
    grr = reaction_gene_association(model, original_rid)[iso_num]
    pstoich = model.enzymedata.reaction_protein_stoichiometry[original_rid][iso_num]
    kcat =
        contains(rid, "§FOR") ? model.enzymedata.reaction_kcats[original_rid][iso_num][1] :
        model.enzymedata.reaction_kcats[original_rid][iso_num][2]
    for (idx, pid) in enumerate(grr)
        push!(E_components.row_idxs, first(indexin([pid], protein_ids)))
        push!(E_components.col_idxs, col_idx)
        push!(E_components.coeffs, -pstoich[idx] / kcat)
    end
end

"""
    _get_proteins_with_kcats(model::GeckoModel)

Return all protein (gene ids) that have a kcat from `model` based on `reaction_kcats` field.
Assume that if a reaction has a kcat then each isozyme has a kcat.
"""
function _get_proteins_with_kcats(model::GeckoModel)
    unique(
        vcat(
            vcat(
                [
                    reaction_gene_association(model.smodel, rid) for
                    rid in reactions(model.smodel) if haskey(model.enzymedata.reaction_kcats, rid)
                ]...,
            )...,
        ),
    )
end
