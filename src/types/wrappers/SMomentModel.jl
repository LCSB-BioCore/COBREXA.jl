
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

Use [`make_simplified_enzyme_constrained_model`](@ref) or [`with_simplified_enzyme_constrained`](@ref) to construct the
models.

The model is constructed as follows:
- stoichiometry of the original model is retained as much as possible, but
  enzymatic reations are split into forward and reverse parts (marked by a
  suffix like `...#forward` and `...#reverse`),
- coupling is added to simulate a virtual metabolite "enzyme capacity", which
  is consumed by all enzymatic reactions at a rate given by enzyme mass divided
  by the corresponding kcat,
- the total consumption of the enzyme capacity is constrained to a fixed
  maximum.

The `SimplifiedEnzymeConstrainedModel` structure contains a worked-out representation of the
optimization problem atop a wrapped [`AbstractMetabolicModel`](@ref), in particular the
separation of certain reactions into unidirectional forward and reverse parts,
an "enzyme capacity" required for each reaction, and the value of the maximum
capacity constraint. Original coupling in the inner model is retained.

In the structure, the field `columns` describes the correspondence of stoichiometry
columns to the stoichiometry and data of the internal wrapped model, and
`total_enzyme_capacity` is the total bound on the enzyme capacity consumption
as specified in sMOMENT algorithm.

This implementation allows easy access to fluxes from the split reactions
(available in `variables(model)`), while the original "simple" reactions from
the wrapped model are retained as [`reactions`](@ref). All additional constraints
are implemented using [`coupling`](@ref) and [`coupling_bounds`](@ref).

# Fields
$(TYPEDFIELDS)
"""
struct SimplifiedEnzymeConstrainedModel <: AbstractModelWrapper
    columns::Vector{_SimplifiedEnzymeConstrainedColumn}
    total_enzyme_capacity::Float64

    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(model::SimplifiedEnzymeConstrainedModel) = model.inner

"""
$(TYPEDSIGNATURES)

Return a stoichiometry of the [`SimplifiedEnzymeConstrainedModel`](@ref). The enzymatic reactions
are split into unidirectional forward and reverse ones.
"""
Accessors.stoichiometry(model::SimplifiedEnzymeConstrainedModel) =
    stoichiometry(model.inner) * simplified_enzyme_constrained_column_reactions(model)

"""
$(TYPEDSIGNATURES)

Reconstruct an objective of the [`SimplifiedEnzymeConstrainedModel`](@ref).
"""
Accessors.objective(model::SimplifiedEnzymeConstrainedModel) =
    simplified_enzyme_constrained_column_reactions(model)' * objective(model.inner)

"""
$(TYPEDSIGNATURES)

Returns the internal reactions in a [`SimplifiedEnzymeConstrainedModel`](@ref) (these may be split
to forward- and reverse-only parts; reactions IDs are mangled accordingly with
suffixes).
"""
Accessors.variables(model::SimplifiedEnzymeConstrainedModel) =
    let inner_reactions = variables(model.inner)
        [
            simplified_enzyme_constrained_reaction_name(
                inner_reactions[col.reaction_idx],
                col.direction,
            ) for col in model.columns
        ]
    end

"""
$(TYPEDSIGNATURES)

The number of reactions (including split ones) in [`SimplifiedEnzymeConstrainedModel`](@ref).
"""
Accessors.n_variables(model::SimplifiedEnzymeConstrainedModel) = length(model.columns)

"""
$(TYPEDSIGNATURES)

Return the variable bounds for [`SimplifiedEnzymeConstrainedModel`](@ref).
"""
Accessors.bounds(model::SimplifiedEnzymeConstrainedModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in [`SimplifiedEnzymeConstrainedModel`](@ref) to the original
fluxes in the wrapped model.
"""
Accessors.reaction_variables(model::SimplifiedEnzymeConstrainedModel) =
    simplified_enzyme_constrained_column_reactions(model)' * reaction_variables(model.inner)

"""
$(TYPEDSIGNATURES)

Return the coupling of [`SimplifiedEnzymeConstrainedModel`](@ref). That combines the coupling of
the wrapped model, coupling for split reactions, and the coupling for the total
enzyme capacity.
"""
Accessors.coupling(model::SimplifiedEnzymeConstrainedModel) = vcat(
    coupling(model.inner) * simplified_enzyme_constrained_column_reactions(model),
    [col.capacity_required for col in model.columns]',
)

"""
$(TYPEDSIGNATURES)

Count the coupling constraints in [`SimplifiedEnzymeConstrainedModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
Accessors.n_coupling_constraints(model::SimplifiedEnzymeConstrainedModel) =
    n_coupling_constraints(model.inner) + 1

"""
$(TYPEDSIGNATURES)

The coupling bounds for [`SimplifiedEnzymeConstrainedModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
Accessors.coupling_bounds(model::SimplifiedEnzymeConstrainedModel) =
    let (iclb, icub) = coupling_bounds(model.inner)
        (vcat(iclb, [0.0]), vcat(icub, [model.total_enzyme_capacity]))
    end
