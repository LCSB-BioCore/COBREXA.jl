"""
    gecko(   
        model::StandardModel,
        optimizer;
        objective_id = "",
        protein_stoichiometry = Dict(),
        protein_masses = Dict(),
        reaction_kcats = Dict(),
        lb_protein_measurements = Dict(),
        ub_protein_measurements = Dict(),
        lb_flux_measurements = Dict(),
        ub_flux_measurements = Dict(),
        total_protein_mass = 0.0,
    )

Perform flux balance analysis on `model` with `optimizer`, using GECKO to
incorporate enzyme capacity and kinetic constraints. See `Sánchez, Benjamín J.,
et al. "Improving the phenotype predictions of a yeast genome‐scale metabolic
model by incorporating enzymatic constraints." Molecular systems biology, 2017.`
for implementation details.

Total enzyme capacity (sum of all enzyme concentrations multiplied by their
molar mass) is constrained by `total_protein_mass`, a unitless mass fraction of
enzyme mass to cell dry mass. The reaction fluxes and protein concentrations can
be bounded by `lb_flux_measurements`, `ub_flux_measurements`,
`lb_protein_measurements`, and `ub_protein_measurements` respectively. Both
lower and upper bounds need to be supplied if a reaction flux is to be bounded,
likewise with protein concentration bounds. The reaction to be optimized is
specified by `objective_id`. Note, since the model uses irreversible reactions
internally, you should append `"§FOR"` for the forward direction and `"§REV"`
for the reverse direction in which ever reaction you want to optimize; this is
    not necesarry for the bound constraints. To optimize anything else, use the
lower level [`gecko_opt_problem`](@ref). Futhermore, `"§"` is reserved for
internal use as a delimiter, no reaction id should contain that character. 
    
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

The function returns a dictionary mapping reaction ids to their fluxes, as well
as a dictionary mapping gene ids to their concentrations. Note, the units depend
on those used in `reaction_kcats` and `protein_masses`. Only the protein and
reaction flux bounds are optional kwargs, all other kwargs must be supplied.
Only reactions with kcats will have enzyme bounds associated with them, but all
isozymes are assumed to have data if data is supplied.

Currently only `modifications` that change attributes of the `optimizer` are 
supported.
"""
function gecko(
    model::StandardModel,
    optimizer;
    objective_id = "",
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    lb_protein_measurements = Dict(),
    ub_protein_measurements = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    total_protein_mass = 0.0,
    sense = MOI.MAX_SENSE,
    modifications = [],
)

    _, E, d, M, h, reaction_map, _, protein_ids = gecko_opt_problem(
        model;
        protein_stoichiometry,
        protein_masses,
        reaction_kcats,
        lb_protein_measurements,
        ub_protein_measurements,
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

    _map_irrev_to_rev_ids(reaction_map, value.(x); protein_ids)
end

"""
    gecko_opt_problem(
        model::StandardModel;
        protein_stoichiometry = Dict(),
        protein_masses = Dict(),
        reaction_kcats = Dict(),
        lb_protein_measurements = Dict(),
        ub_protein_measurements = Dict(),
        lb_flux_measurements = Dict(),
        ub_flux_measurements = Dict(),
        total_protein_mass = 0.0,
    )

Lower level function that returns the matrix form of a model with enzyme capacity 
constraints, in GECKO format, see [`gecko`](@ref) for the higher level function.
```
max/min cᵀ * x
s.t.    E * x = d 
        M * x ≤ h
```
Returns `c, E, d, M, h, reaction_map, metabolite_map, protein_ids`, where 
`reaction_map` shows the order of the columns (reactions) in `E`. Proteins 
are ordered according to `protein_ids`, and follow after reactions. Use 
[`_map_irrev_to_rev_ids`](@ref) to map the solution of an optimization 
problem back to the original model's name space.

Format of arguments are always in order of grr for each reaction `rxn_id`:
1) protein_stoichiometry: `Dict(rxn_id => [[1,2,1,1]])` 
2) protein_masses: `Dict(p_id => [mm, ...])` in units of kDa
3) reaction_kcat: `Dict(rxn_id => [[kcat_for, kcat_rev]])` for each complex
    

Assumptions:
1) Each isozyme has a kcat (forward and reverse) for each reaction it catalyzes
2) Only reactions with kcats have enzyme constraints
3) Both `lb_flux_measurements` and `ub_flux_measurements` have the same keys

Notes:
1) The objective vector, `c` is not set
2) The parameters are the kcats and the total protein measurement
3) The symbol `§` is a reserved delimiter, do not use it in reaction or metabolite ids
"""
function gecko_opt_problem(
    model::StandardModel;
    protein_stoichiometry = Dict(),
    protein_masses = Dict(),
    reaction_kcats = Dict(),
    lb_protein_measurements = Dict(),
    ub_protein_measurements = Dict(),
    lb_flux_measurements = Dict(),
    ub_flux_measurements = Dict(),
    total_protein_mass = 0.0,
)
    S, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        _build_irreversible_stoichiometric_matrix(model)

    #: find all gene products that have kcats associated with them
    protein_ids = _get_proteins_with_kcats(model, reaction_kcats)

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
        !haskey(reaction_kcats, original_rid) && continue

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
            protein_stoichiometry,
            reaction_kcats,
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

    #: need to set objective reaction outside
    c = spzeros(n_vars)

    #: inequality constraints
    M, h = _gecko_build_inequality_constraints(
        lb_protein_measurements,
        ub_protein_measurements,
        protein_ids,
        protein_masses,
        n_reactions,
        n_proteins,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
        total_protein_mass,
    )

    return c, E, d, M, h, reaction_map, metabolite_map, protein_ids
end

"""
    _gecko_build_inequality_constraints(
        lb_protein_measurements,
        ub_protein_measurements,
        protein_ids,
        protein_masses,
        n_reactions,
        n_proteins,
        lb_flux_measurements,
        ub_flux_measurements,
        lb_fluxes,
        ub_fluxes,
        reaction_map,
        total_protein_mass,
    )

Helper function to build inequality constraints. Returns the inequality constraint in matrix format.
"""
function _gecko_build_inequality_constraints(
    lb_protein_measurements,
    ub_protein_measurements,
    protein_ids,
    protein_masses,
    n_reactions,
    n_proteins,
    lb_flux_measurements,
    ub_flux_measurements,
    lb_fluxes,
    ub_fluxes,
    reaction_map,
    total_protein_mass,
)
    #: inequality lhs
    mw_proteins = [protein_masses[pid] for pid in protein_ids]
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
    for original_rid in keys(lb_flux_measurements) # only constrain if measurement available
        lb = lb_flux_measurements[original_rid]
        ub = ub_flux_measurements[original_rid]
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
        haskey(lb_protein_measurements, pid) ? lb_protein_measurements[pid] : 0.0 for
        pid in protein_ids
    ]
    ub_proteins = [
        haskey(ub_protein_measurements, pid) ? ub_protein_measurements[pid] : 10_000.0
        for pid in protein_ids
    ]

    h = Array([-lb_fluxes; ub_fluxes; -lb_proteins; ub_proteins; total_protein_mass])

    return M, h
end
