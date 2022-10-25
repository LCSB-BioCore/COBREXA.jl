
"""
$(TYPEDEF)

A helper type that describes the contents of [`SMomentModel`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
struct _SMomentColumn
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

Use [`make_smoment_model`](@ref) or [`with_smoment`](@ref) to construct the
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

The `SMomentModel` structure contains a worked-out representation of the
optimization problem atop a wrapped [`AbstractMetabolicModel`](@ref), in particular the
separation of certain reactions into unidirectional forward and reverse parts,
an "enzyme capacity" required for each reaction, and the value of the maximum
capacity constraint. Original coupling in the inner model is retained.

In the structure, the field `columns` describes the correspondence of stoichiometry
columns to the stoichiometry and data of the internal wrapped model, and
`total_enzyme_capacity` is the total bound on the enzyme capacity consumption
as specified in sMOMENT algorithm.

This implementation allows easy access to fluxes from the split reactions
(available in `reactions(model)`), while the original "simple" reactions from
the wrapped model are retained as [`fluxes`](@ref). All additional constraints
are implemented using [`coupling`](@ref) and [`coupling_bounds`](@ref).

# Fields
$(TYPEDFIELDS)
"""
struct SMomentModel <: ModelWrapper
    columns::Vector{_SMomentColumn}
    total_enzyme_capacity::Float64

    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(model::SMomentModel) = model.inner

"""
$(TYPEDSIGNATURES)

Return a stoichiometry of the [`SMomentModel`](@ref). The enzymatic reactions
are split into unidirectional forward and reverse ones.
"""
Accessors.stoichiometry(model::SMomentModel) =
    stoichiometry(model.inner) * smoment_column_reactions(model)

"""
$(TYPEDSIGNATURES)

Reconstruct an objective of the [`SMomentModel`](@ref).
"""
Accessors.objective(model::SMomentModel) =
    smoment_column_reactions(model)' * objective(model.inner)

"""
$(TYPEDSIGNATURES)

Returns the internal reactions in a [`SMomentModel`](@ref) (these may be split
to forward- and reverse-only parts; reactions IDs are mangled accordingly with
suffixes).
"""
Accessors.reactions(model::SMomentModel) =
    let inner_reactions = reactions(model.inner)
        [
            smoment_reaction_name(inner_reactions[col.reaction_idx], col.direction) for
            col in model.columns
        ]
    end

"""
$(TYPEDSIGNATURES)

The number of reactions (including split ones) in [`SMomentModel`](@ref).
"""
Accessors.n_reactions(model::SMomentModel) = length(model.columns)

"""
$(TYPEDSIGNATURES)

Return the variable bounds for [`SMomentModel`](@ref).
"""
Accessors.bounds(model::SMomentModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in [`SMomentModel`](@ref) to the original
fluxes in the wrapped model.
"""
Accessors.reaction_flux(model::SMomentModel) =
    smoment_column_reactions(model)' * reaction_flux(model.inner)

"""
$(TYPEDSIGNATURES)

Return the coupling of [`SMomentModel`](@ref). That combines the coupling of
the wrapped model, coupling for split reactions, and the coupling for the total
enzyme capacity.
"""
Accessors.coupling(model::SMomentModel) = vcat(
    coupling(model.inner) * smoment_column_reactions(model),
    [col.capacity_required for col in model.columns]',
)

"""
$(TYPEDSIGNATURES)

Count the coupling constraints in [`SMomentModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
Accessors.n_coupling_constraints(model::SMomentModel) =
    n_coupling_constraints(model.inner) + 1

"""
$(TYPEDSIGNATURES)

The coupling bounds for [`SMomentModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
Accessors.coupling_bounds(model::SMomentModel) =
    let (iclb, icub) = coupling_bounds(model.inner)
        (vcat(iclb, [0.0]), vcat(icub, [model.total_enzyme_capacity]))
    end
