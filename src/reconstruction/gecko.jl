

"""
    with_gecko(; kwargs...)

Specifies a model variant which adds extra semantics of the Gecko algorithm,
giving a [`GeckoModel`](@ref). The arguments are forwarded to
[`make_gecko_model`](@ref). Intended for usage with [`screen`](@ref).
"""
with_gecko(; kwargs...) = model -> make_gecko_model(model; kwargs...)

"""
    change_bound!(model::GeckoModel, id; lower=nothing, upper=nothing)

Change the bound of variable in `model`. Does not change the bound if respective
bound is `nothing`. Note, for `GeckoModel`s, if the model used to construct the
`GeckoModel` has irreversible reactions, then these reactions will be
permanently irreversible in the model, i.e. changing their bounds to make them
reversible will have no effect.
"""
function change_bound!(model::GeckoModel, id; lower = nothing, upper = nothing)
    #TODO
end

"""
    change_bounds!(model::GeckoModel, ids; lower=fill(nothing, length(ids)), upper=fill(nothing, length(ids)))

Change the bounds of multiple variables in `model` simultaneously. See
[`change_bound`](@ref) for details.
"""
function change_bounds!(
    model::GeckoModel,
    ids;
    lower = fill(nothing, length(ids)),
    upper = fill(nothing, length(ids)),
)
    #TODO
end
