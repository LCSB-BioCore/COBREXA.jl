
"""
$(TYPEDSIGNATURES)

Run a basic enzyme constrained flux balance analysis on `model`. The enzyme
model is parameterized by `reaction_isozymes`, which is a mapping of reaction
IDs (those used in the fluxes of the model) to named [`Isozyme`](@ref)s.
Additionally, `gene_molar_masses` and `capacity_limitations` should be supplied.
The latter is a vector of tuples, where each tuple represents a distinct bound
as `(bound_id, genes_in_bound, protein_mass_bound)`. Typically, `model` has
bounded exchange reactions, which are unnecessary in enzyme constrained models.
Unbound these reactions by listing their IDs in `unconstrain_reactions`, which
makes them reversible. Optimization `modifications` are directly forwarded to
[`optimized_constraints`](@ref).

In the event that your model requires more complex build steps, consider
constructing it manually by using [`add_enzyme_constraints!`](@ref).

Returns a tree with the optimization solution or `nothing` if the model cannot
be solved.
"""
function enzyme_constrained_flux_balance_analysis(
    model::A.AbstractFBCModel,
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_molar_masses::Dict{String,Float64},
    capacity_limitations::Vector{Tuple{String,Vector{String},Float64}};
    optimizer,
    unconstrain_reactions = String[],
    modifications = [],
)
    m = build_flux_balance_model(model)

    # create enzyme variables
    m += :enzymes^enzyme_variables(model)

    m = add_enzyme_constraints!(
        m,
        reaction_isozymes,
        gene_molar_masses,
        capacity_limitations,
    )

    for rid in Symbol.(unconstrain_reactions)
        m.fluxes[rid].bound = (-1000.0, 1000.0)
    end

    optimized_constraints(m; objective = m.objective.value, optimizer, modifications)
end

export enzyme_constrained_flux_balance_analysis
