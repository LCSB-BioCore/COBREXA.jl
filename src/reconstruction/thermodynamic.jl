"""
$(TYPEDSIGNATURES)

Construct a [`MaxMinDrivingForceModel`](@ref) so that max min driving force
analysis can be performed on `model`.
"""
make_max_min_driving_force_model(model::AbstractMetabolicModel; kwargs...) =
    MaxMinDrivingForceModel(; inner = model, kwargs...)
