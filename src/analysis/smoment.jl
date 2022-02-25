"""
    smoment(
        model::StandardModel,
        optimizer;
        objective_id = "",
        protein_stoichiometry = Dict(),
        protein_masses = Dict(),
        reaction_kcats = Dict(),
        lb_flux_measurements = Dict(),
        ub_flux_measurements = Dict(),
        total_protein_mass = 0.0,
        sense = MOI.MAX_SENSE,
        modifications = [],
    )

Perform enzyme capacity constrained flux balance analysis on `model` with
`optimizer` using the SMOMENT algorithm, see `Bekiaris, Pavlos Stephanos, and
Steffen Klamt. "Automatic construction of metabolic models with enzyme
constraints." BMC bioinformatics, 2020.` for implementation details. 

SMOMENT is a direct simplification of GECKO (despite it being named after the
MOMENT algorithm). Total enzyme capacity (sum of all enzyme concentrations
multiplied by their molar mass) is constrained by `total_protein_mass`, a
unitless mass fraction of enzyme mass to cell dry mass. The reaction fluxes can
be bounded by `lb_flux_measurements`, `ub_flux_measurements`. Both lower and
upper bounds need to be supplied if a reaction flux is to be bounded. The
reaction to be optimized is specified by `objective_id`. Note, since the
model uses irreversible reactions internally, you should append `"§FOR"` for the
forward direction and `"§REV"` for the reverse direction in which ever reaction
you want to optimize; this is not necesarry for the bound constraints. To
optimize anything else, use the lower level [`smoment_opt_problem`](@ref).
Futhermore, `"§"` is reserved for internal use as a delimiter, no reaction id
should contain that character. Also note, SMOMENT assumes that each reaction only has 
a single enzyme (one GRR) associated with it. It is required that a model be modified to
ensure that this condition is met. For ease-of-use, [`remove_slow_isozymes!`](@ref) is 
supplied to effect this. 
    
The protein masses (in molar mass units) for each gene in the model should also
be supplied through `protein_masses`. The format is a dictionary of gene ids
mapped to molar masses. Additionally, the reaction turnover numbers (catalytic
constants, kcats) are supplied through `reaction_kcats`, which is also a
dictionary mapping reaction ids to kcats of each isozyme encoded by the
reaction's gene reaction rule. Each isozyme should have both a forward and
reverse kcat, so `reaction_kcats = Dict(rid => [[k1f, k1r], [k2f, k2r]], ...)`
for `rid` with two isozymes. Finally, the stoichiometry of each isozyme needs to
be supplied by `protein_stoichiometry`. The format is also a dictionary mapping
gene ids returned by [`reaction_gene_association`](@ref) to their stoichiometry,
e.g. `protein_stoichiometry = Dict(rid => [[1,1],[1,2]],...)` implies that the
first isozyme of `rid` is composed of two subunits, each present once in the
protein, while the second isozyme is composed of two subunits, but the second
subunit is present twice in the isozyme.

The function returns a dictionary mapping reaction ids to their fluxes. Note,
the units depend on those used in `reaction_kcats` and `protein_masses`. Only
the protein and reaction flux bounds are optional kwargs, all other kwargs must
be supplied. Only reactions with kcats will have enzyme bounds associated with
them, but all isozymes are assumed to have data if data is supplied.

Currently only `modifications` that change attributes of the `optimizer` are 
supported. 
"""
function smoment(
    model::StandardModel,
    optimizer;
    objective_id = "",
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    total_protein_mass = 0.0,
    sense = MOI.MAX_SENSE,
    modifications = [],
)

    _, E, d, M, h, reaction_map, _ = smoment_opt_problem(
        model;
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        lb_flux_measurements,
        ub_flux_measurements,
        total_protein_mass,
    )

    opt_model = Model(optimizer)
    x = @variable(opt_model, x[1:size(E, 2)])
    bid = reaction_map[objective_id]
    @objective(opt_model, sense, x[bid])
    @constraint(opt_model, E * x .== d)
    @constraint(opt_model, M * x .<= h)

    # apply the modifications, if any
    for mod in modifications
        mod(nothing, opt_model)
    end

    optimize!(opt_model)

    _map_irrev_to_rev_ids(reaction_map, value.(x))

end

"""
    smoment_opt_problem(
        model::StandardModel;
        protein_stoichiometry = Dict(),
        protein_masses = Dict(),
        reaction_kcats = Dict(),
        lb_flux_measurements = Dict(),
        ub_flux_measurements = Dict(),
        total_protein_mass = 0.0,
    )

Lower level function that returns the matrix form of a model with enzyme capacity 
constraints, in SMOMENT format, see [`smoment`](@ref) for the higher level function.
```
max/min cᵀ * x
s.t.    E * x = d 
        M * x ≤ h
```
Returns `c, E, d, M, h, reaction_map, metabolite_map`, where `reaction_map`
shows the order of the columns (reactions) in `E`. Use
[`_map_irrev_to_rev_ids`](@ref) to map the solution of an optimization problem
back to the original model's name space. Note, this function implements the most
basic version of SMOMENT, i.e. you cannot limit the concentration of any protein
(use [`gecko`](@ref) for that). Importantly, this function assumes that a
preprocessing step has been performed that changes the model so that each
reaction only has one GRR corresponding to the fastest isozyme. For this
preprocessing step, use [`remove_slow_isozymes!`](@ref).

Format of arguments are always in order of grr for each reaction `rxn_id`:
1) protein_stoichiometry: `Dict(rxn_id => [[1,2,1,1]])` 
2) protein_masses: `Dict(p_id => [mm, ...])` in units of kDa
3) reaction_kcat: `Dict(rxn_id => [[kcat_for, kcat_rev]])` NOTE: no isozymes.

Assumptions:
1) No isozymes.
2) Both `lb_flux_measurements` and `ub_flux_measurements` have the same keys

Notes:
1) The objective vector, `c` is not set
2) The parameters are the kcats and the total protein measurement
3) The symbol `§` is a reserved delimiter, do not use it in reaction or metabolite ids
"""
function smoment_opt_problem(
    model::StandardModel;
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    total_protein_mass = 0.0,
)

    S, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        _build_irreversible_stoichiometric_matrix(model)

    #: size of resultant model
    n_reactions = size(S, 2)
    n_metabolites = size(S, 1)
    n_vars = n_reactions + 1

    #: equality lhs
    Se = zeros(1, n_reactions)

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        !haskey(reaction_kcats, original_rid) && continue
        # these entries have kcats, only one GRR by assumption
        grr = first(reaction_gene_association(model, original_rid))
        pstoich = first(protein_stoichiometry[original_rid])
        mw = dot(pstoich, [protein_masses[gid] for gid in grr])
        kcat =
            contains(rid, "§FOR") ? first(reaction_kcats[original_rid])[1] :
            first(reaction_kcats[original_rid])[2]
        Se[1, col_idx] = -mw / kcat
    end

    E = [
        S zeros(n_metabolites, 1)
        Se 1.0
    ]

    # #: equality rhs
    d = zeros(n_metabolites + 1)

    # #: need to set objective reaction outside
    c = spzeros(n_vars)

    #: inequality constraints
    M, h = _smoment_build_inequality_constraints(
        n_reactions,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
        total_protein_mass,
    )

    return c, E, d, M, h, reaction_map, metabolite_map
end


"""
    _smoment_build_inequality_constraints(
        n_reactions,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
    )

Helper function to return functions describing the inequality 
constraints for smoment.
"""
function _smoment_build_inequality_constraints(
    n_reactions,
    lb_flux_measurements,
    ub_flux_measurements,
    lb_fluxes,
    ub_fluxes,
    reaction_map,
    total_protein_mass,
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
    for original_rid in keys(lb_flux_measurements) # only constrain if measurement available
        lb = lb_flux_measurements[original_rid]
        ub = ub_flux_measurements[original_rid]
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

    h = Array([-lb_fluxes; ub_fluxes; total_protein_mass])

    return M, h
end
