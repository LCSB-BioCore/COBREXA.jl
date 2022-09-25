
"""
$(TYPEDEF)

A helper type that describes the contents of [`SMomentModel`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
struct _smoment_column
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

# Fields
$(TYPEDFIELDS)
"""
struct SMomentModel <: ModelWrapper
    columns::Vector{_smoment_column}
    total_enzyme_capacity::Float64

    inner::MetabolicModel
end

# smoment 

"""
$(TYPEDSIGNATURES)

Internal helper for systematically naming reactions in [`SMomentModel`](@ref).
"""
_smoment_reaction_name(original_name::String, direction::Int) =
    direction == 0 ? original_name :
    direction > 0 ? "$original_name#forward" : "$original_name#reverse"

"""
$(TYPEDSIGNATURES)

Retrieve a utility mapping between reactions and split reactions; rows
correspond to "original" reactions, columns correspond to "split" reactions.
"""
_smoment_column_reactions(model::SMomentModel) = sparse(
    [col.reaction_idx for col in model.columns],
    1:length(model.columns),
    [col.direction >= 0 ? 1 : -1 for col in model.columns],
    n_reactions(model.inner),
    length(model.columns),
)