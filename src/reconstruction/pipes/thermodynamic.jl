"""
$(TYPEDSIGNATURES)

A pipe-able function that specifies a model variant that solves the max-min
driving force problem. Calls [`make_max_min_driving_force_model`](@ref)
internally.
"""
with_max_min_driving_force_analysis(args...; kwargs...) =
    m -> make_max_min_driving_force_model(m, args...; kwargs...)
