
"""
$(TYPEDEF)

A helper type that describes the contents of [`SimplifiedEnzymeConstrainedModel`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
struct _SimplifiedEnzymeConstrainedColumn
    reaction_idx::Int # number of the corresponding reaction in the inner model
    direction::Int # 0 if "as is" and unique, -1 if reverse-only part, 1 if forward-only part
    lb::Float64 # must be 0 if the reaction is unidirectional (if direction!=0)
    ub::Float64
    capacity_required::Float64 # must be 0 for bidirectional reactions (if direction==0)
end

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
- coupling is added to simulate a virtual metabolite "enzyme capacity", which is
  consumed by all enzymatic reactions at a rate given by enzyme mass divided by
  the corresponding kcat,
- the total consumption of the enzyme capacity is constrained to a fixed
  maximum.

The `SimplifiedEnzymeConstrainedModel` structure contains a worked-out
representation of the optimization problem atop a wrapped
[`AbstractMetabolicModel`](@ref), in particular the separation of certain
reactions into unidirectional forward and reverse parts (which changes the
stoichiometric matrix), an "enzyme capacity" required for each reaction, and the
value of the maximum capacity constraint. Original coupling in the inner model
is retained.

In the structure, the field `columns` describes the correspondence of
stoichiometry columns to the stoichiometry and data of the internal wrapped
model, and `total_gene_product_mass_bound` is the total bound on the enzyme
capacity consumption as specified in sMOMENT algorithm.

This implementation allows easy access to fluxes from the split reactions
(available in `variables(model)`), while the original "simple" reactions from
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
    columns::Vector{_SimplifiedEnzymeConstrainedColumn}
    total_gene_product_mass_bound::Float64

    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(model::SimplifiedEnzymeConstrainedModel) = model.inner

Accessors.stoichiometry(model::SimplifiedEnzymeConstrainedModel) =
    stoichiometry(model.inner) * simplified_enzyme_constrained_column_reactions(model)

Accessors.objective(model::SimplifiedEnzymeConstrainedModel) =
    simplified_enzyme_constrained_column_reactions(model)' * objective(model.inner)

Accessors.variables(model::SimplifiedEnzymeConstrainedModel) =
    let inner_reactions = variables(model.inner)
        [
            simplified_enzyme_constrained_reaction_name(
                inner_reactions[col.reaction_idx],
                col.direction,
            ) for col in model.columns
        ]
    end

Accessors.n_variables(model::SimplifiedEnzymeConstrainedModel) = length(model.columns)

Accessors.bounds(model::SimplifiedEnzymeConstrainedModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in
[`SimplifiedEnzymeConstrainedModel`](@ref) to the original fluxes in the
wrapped model (as a matrix).
"""
Accessors.reaction_variables_matrix(model::SimplifiedEnzymeConstrainedModel) =
    simplified_enzyme_constrained_column_reactions(model)' *
    reaction_variables_matrix(model.inner)

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in
[`SimplifiedEnzymeConstrainedModel`](@ref) to the original fluxes in the
wrapped model.
"""
Accessors.reaction_variables(model::SimplifiedEnzymeConstrainedModel) =
    Accessors.Internal.make_mapping_dict(
        variables(model),
        semantics(model),
        reaction_variables_matrix(model),
    ) # TODO currently inefficient

Accessors.coupling(model::SimplifiedEnzymeConstrainedModel) = vcat(
    coupling(model.inner) * simplified_enzyme_constrained_column_reactions(model),
    [col.capacity_required for col in model.columns]',
)

Accessors.n_coupling_constraints(model::SimplifiedEnzymeConstrainedModel) =
    n_coupling_constraints(model.inner) + 1

Accessors.coupling_bounds(model::SimplifiedEnzymeConstrainedModel) =
    let (iclb, icub) = coupling_bounds(model.inner)
        (vcat(iclb, [0.0]), vcat(icub, [model.total_gene_product_mass_bound]))
    end
