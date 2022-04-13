

"""
    with_smoment(; kwargs...)

Specifies a model variant which adds extra semantics of the sMOMENT algorithm,
giving a [`SMomentModel`](@ref). The arguments are forwarded to
[`make_smoment_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_smoment(; kwargs...) = model -> make_smoment_model(model; kwargs...)

"""
    change_bound!(model::SMomentModel, id; lower=nothing, upper=nothing)

Change the bound of variable in `model`. Does not change the bound if respective
bound is `nothing`. Note, for `SMomentModel`s, if the model used to construct the
`SMomentModel` has irreversible reactions, then these reactions will be
permanently irreversible in the model, i.e. changing their bounds to make them
reversible will have no effect.
"""
function change_bound!(model::SMomentModel, id; lower = nothing, upper = nothing)
    # TODO
end

"""
    change_bounds!(model::SMomentModel, ids; lower=fill(nothing, length(ids)), upper=fill(nothing, length(ids)))

Change the bounds of multiple variables in `model` simultaneously. See
[`change_bound`](@ref) for details.
"""
function change_bounds!(
    model::SMomentModel,
    ids;
    lower = fill(nothing, length(ids)),
    upper = fill(nothing, length(ids)),
)
    # TODO
end
