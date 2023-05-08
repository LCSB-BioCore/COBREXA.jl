"""
$(TYPEDEF)

An enzyme-capacity-constrained model using sMOMENT algorithm, as described by
*Bekiaris, Pavlos Stephanos, and Steffen Klamt, "Automatic construction of
metabolic models with enzyme constraints" BMC bioinformatics, 2020*.

Use [`make_simplified_enzyme_constrained_model`](@ref) or
[`with_simplified_enzyme_constraints`](@ref) to construct the models.

The model is constructed as follows:
- stoichiometry of the original model is retained as much as possible, but
  enzymatic reations are split into forward and reverse parts (marked by a
  suffix like `...#forward` and `...#reverse`),
- coupling is added to simulate lumped virtual metabolites that act like enzyme
  capacities. These are consumed by enzymatic reactions at a rate given by
  enzyme mass divided by the corresponding kcat,
- the total consumption of the enzyme capacity bounds is constrained to be less
  than some fixed values.

The `SimplifiedEnzymeConstrainedModel` structure contains a worked-out
representation of the optimization problem atop a wrapped
[`AbstractMetabolicModel`](@ref). The internal representation of the model
splits reactions into unidirectional forward and reverse parts (which changes
the stoichiometric matrix).

In the structure, the field `columns` describes the correspondence of
stoichiometry columns to the stoichiometry, and data of the internal wrapped
model. Multiple capacity bounds may be added through
`total_reaction_mass_bounds`. These bounds are connected to the model through
`columns`. Since this algorithm is reaction centered, no enzymes directly appear
in the formulation.

This implementation allows easy access to fluxes from the split reactions
(available in `variable_ids(model)`), while the original "simple" reactions from
the wrapped model are retained as [`reactions`](@ref). All additional
constraints are implemented using [`coupling`](@ref) and
[`coupling_bounds`](@ref).

To implement this wrapper for a model, the accessors [`reaction_isozymes`](@ref)
and [`gene_product_molar_mass`](@ref), need to be available. Additionally, the
model needs to associate [`Isozyme`](@ref)s with reactions. Reactions without
enzymes, or those that should be ignored need to return `nothing` when
[`reaction_isozymes`](@ref) is called on them.

# Fields
$(TYPEDFIELDS)
"""
struct SimplifiedEnzymeConstrainedModel <: AbstractModelWrapper
    columns::Vector{SimplifiedEnzymeConstrainedColumn}
    total_reaction_mass_bounds::Vector{Float64}

    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(model::SimplifiedEnzymeConstrainedModel) = model.inner

Accessors.stoichiometry(model::SimplifiedEnzymeConstrainedModel) =
    stoichiometry(model.inner) * simplified_enzyme_constrained_column_reactions(model)

Accessors.objective(model::SimplifiedEnzymeConstrainedModel) =
    simplified_enzyme_constrained_column_reactions(model)' * objective(model.inner)

Accessors.variable_ids(model::SimplifiedEnzymeConstrainedModel) =
    let inner_reactions = variable_ids(model.inner)
        [
            simplified_enzyme_constrained_reaction_name(
                inner_reactions[col.reaction_idx],
                col.direction,
            ) for col in model.columns
        ]
    end

Accessors.variable_count(model::SimplifiedEnzymeConstrainedModel) = length(model.columns)

Accessors.variable_bounds(model::SimplifiedEnzymeConstrainedModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in
[`SimplifiedEnzymeConstrainedModel`](@ref) to the original fluxes in the
wrapped model (as a matrix).
"""
Accessors.reaction_variables_matrix(model::SimplifiedEnzymeConstrainedModel) =
    reaction_variables_matrix(model.inner) *
    simplified_enzyme_constrained_column_reactions(model) #TODO check

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in
[`SimplifiedEnzymeConstrainedModel`](@ref) to the original fluxes in the
wrapped model.
"""
Accessors.reaction_variables(model::SimplifiedEnzymeConstrainedModel) =
    Accessors.Internal.make_mapping_dict(
        reaction_ids(model.inner),
        variable_ids(model),
        reaction_variables_matrix(model),
    ) # TODO currently inefficient

function Accessors.coupling(model::SimplifiedEnzymeConstrainedModel)
    inner_coupling =
        coupling(model.inner) * simplified_enzyme_constrained_column_reactions(model)

    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (col_idx, col) in enumerate(model.columns)
        for row_idx in col.capacity_bound_idxs
            push!(J, col_idx)
            push!(I, row_idx)
            push!(V, col.capacity_contribution)
        end
    end

    capacity_coupling =
        sparse(I, J, V, length(model.total_reaction_mass_bounds), length(model.columns))

    return [
        inner_coupling
        capacity_coupling
    ]
end

Accessors.n_coupling_constraints(model::SimplifiedEnzymeConstrainedModel) =
    n_coupling_constraints(model.inner) + length(model.total_reaction_mass_bounds)

Accessors.coupling_bounds(model::SimplifiedEnzymeConstrainedModel) =
    let (iclb, icub) = coupling_bounds(model.inner)
        (
            vcat(iclb, zeros(length(model.total_reaction_mass_bounds))),
            vcat(icub, model.total_reaction_mass_bounds),
        )
    end
