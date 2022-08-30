
"""
    struct _smoment_column

A helper type that describes the contents of [`SMomentModel`](@ref)s.
"""
struct _smoment_column
    reaction_idx::Int # number of the corresponding reaction in the inner model
    direction::Int # 0 if "as is" and unique, -1 if reverse-only part, 1 if forward-only part
    lb::Float64 # must be 0 if the reaction is unidirectional (if direction!=0)
    ub::Float64
    capacity_required::Float64 # must be 0 for bidirectional reactions (if direction==0)
end

"""
    struct SMomentModel <: ModelWrapper

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
optimization problem atop a wrapped [`MetabolicModel`](@ref), in particular the
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
"""
@with_repr struct SMomentModel <: ModelWrapper
    columns::Vector{_smoment_column}
    total_enzyme_capacity::Float64

    inner::MetabolicModel
end

unwrap_model(model::SMomentModel) = model.inner

"""
    stoichiometry(model::SMomentModel)

Return a stoichiometry of the [`SMomentModel`](@ref). The enzymatic reactions
are split into unidirectional forward and reverse ones.
"""
stoichiometry(model::SMomentModel) =
    stoichiometry(model.inner) * _smoment_column_reactions(model)

"""
    objective(model::SMomentModel)

Reconstruct an objective of the [`SMomentModel`](@ref).
"""
objective(model::SMomentModel) = _smoment_column_reactions(model)' * objective(model.inner)

"""
    reactions(model::SMomentModel)

Returns the internal reactions in a [`SMomentModel`](@ref) (these may be split
to forward- and reverse-only parts; reactions IDs are mangled accordingly with
suffixes).
"""
reactions(model::SMomentModel) =
    let inner_reactions = reactions(model.inner)
        [
            _smoment_reaction_name(inner_reactions[col.reaction_idx], col.direction) for
            col in model.columns
        ]
    end

"""
    n_reactions(model::SMomentModel)

The number of reactions (including split ones) in [`SMomentModel`](@ref).
"""
n_reactions(model::SMomentModel) = length(model.columns)

"""
    bounds(model::SMomentModel)

Return the variable bounds for [`SMomentModel`](@ref).
"""
bounds(model::SMomentModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
    reaction_flux(model::SMomentModel)

Get the mapping of the reaction rates in [`SMomentModel`](@ref) to the original
fluxes in the wrapped model.
"""
reaction_flux(model::SMomentModel) =
    _smoment_column_reactions(model)' * reaction_flux(model.inner)

"""
    coupling(model::SMomentModel)

Return the coupling of [`SMomentModel`](@ref). That combines the coupling of
the wrapped model, coupling for split reactions, and the coupling for the total
enzyme capacity.
"""
coupling(model::SMomentModel) = vcat(
    coupling(model.inner) * _smoment_column_reactions(model),
    [col.capacity_required for col in model.columns]',
)

"""
    n_coupling_constraints(model::SMomentModel)

Count the coupling constraints in [`SMomentModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
n_coupling_constraints(model::SMomentModel) = n_coupling_constraints(model.inner) + 1

"""
    coupling_bounds(model::SMomentModel)

The coupling bounds for [`SMomentModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
coupling_bounds(model::SMomentModel) =
    let (iclb, icub) = coupling_bounds(model.inner)
        (vcat(iclb, [0.0]), vcat(icub, [model.total_enzyme_capacity]))
    end
