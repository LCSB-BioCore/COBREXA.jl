
"""
    struct _smoment_column

A helper type that describes the contents of [`SMomentModel`](@ref)s.
"""
struct _smoment_column
    reaction_id::Int # number of the corresponding reaction in the inner model
    direction::Int # 0 if "as is" and unique, -1 if reverse-only part, 1 if forward-only part
    coupling_row::Int # number of row in the coupling (0 if direction==0)
    lb::Float64 # must be 0 if the reaction is unidirectional (if direction!=0)
    ub::Float64
    capacity_required::Float64 # must be 0 for bidirectional reactions (if direction==0)
end

# TODO fix the docstring
"""
    struct SMomentModel <: ModelWrapper

Construct an enzyme-capacity constrained model using sMOMENT algorithm, as
described by *Bekiaris, Pavlos Stephanos, and Steffen Klamt, "Automatic
construction of metabolic models with enzyme constraints" BMC bioinformatics,
2020*.

Use [`make_smoment_model`](@ref) or [`with_smoment`](@ref) to construct the
models.

The model is constructed as follows:
- stoichiometry of the original model is retained as much as possible, but
  enzymatic reations are split into forward and reverse parts (marked by a
  suffix like `...#forward` and `...#reverse`),
- sums of forward and reverse reaction pair fluxes are constrained accordingly
  to the original model,
- stoichiometry is expanded by a virtual metabolite "enzyme capacity" which is
  consumed by all enzymatic reactions at a rate given by enzyme mass divided by
  the corresponding kcat,
- the total consumption of the enzyme capacity is constrained by a fixed
  maximum.

The `SMomentModel` structure contains a worked-out representation of the
optimization problem atop a wrapped [`MetabolicModel`](@ref), in particular the
separation of certain reactions into unidirectional forward and reverse parts,
the grouping of these reactions together into virtual "arm" reactions constrained
by bounds from the inner model, an "enzyme capacity" required for each
reaction, and the value of the maximum capacity constraint.

In the structure, field `columns` describes the correspondence of stoichiometry
columns to the stoichiometry and data of the internal wrapped model; field
`coupling_row_reaction` maps the generated coupling constraints to reaction
indexes in the wrapped model, and `total_enzyme_capacity` is the total bound on
the enzyme capacity consumption as specified in sMOMENT algorithm.

This implementation allows easy access to fluxes from the split reactions
(available in `reactions(model)`), while the original "simple" reactions from
the wrapped model are retained as [`fluxes`](@ref). All additional constraints
are implemented using [`coupling`](@ref) and [`coupling_bounds`](@ref).
Original coupling is retained.
"""
struct SMomentModel <: ModelWrapper
    columns::Vector{_smoment_column}
    coupling_row_reaction::Vector{Int}
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
to forward- and reverse-only parts; reactions IDs mangled accordingly with
suffixes).
"""
reactions(model::SMomentModel) =
    let inner_reactions = reactions(model.inner)
        [
            _smoment_reaction_name(inner_reactions[col.reaction_id], col.direction) for
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
    reaction_flux(model.inner) * _smoment_column_reactions(model)

"""
    coupling(model::SMomentModel)

Return the coupling of [`SMomentModel`](@ref). That combines the coupling of
the wrapped model, coupling for split reactions, and the coupling for the total
enzyme capacity.
"""
coupling(model::SMomentModel) = vcat(
    coupling(model.inner) * _smoment_column_reactions(model),
    _smoment_reaction_coupling(model),
    [col.capacity_required for col in model.columns]',
)

"""
    n_coupling_constraints(model::SMomentModel)

Count the coupling constraints in [`SMomentModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
n_coupling_constraints(model::SMomentModel) =
    n_coupling_constraints(model.inner) + _smoment_n_reaction_couplings(model) + 1

"""
    coupling_bounds(model::SMomentModel)

The coupling bounds for [`SMomentModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
function coupling_bounds(model::SMomentModel)
    (ilb, iub) = coupling_bounds(model.inner)
    (rlb, rub) = _smoment_reaction_coupling_bounds(model)
    return (vcat(ilb, rlb, [0.0]), vcat(iub, rub, [model.total_enzyme_capacity]))
end
